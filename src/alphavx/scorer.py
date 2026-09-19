"""Core variant scoring engine wrapping the AlphaGenome SDK.

Handles single-variant scoring, batch processing with retry logic,
and integration with the result cache.
"""

from __future__ import annotations

import logging
import os
import time
from collections.abc import Callable
from typing import TYPE_CHECKING

import pandas as pd

from .config import Config
from .vcf_parser import VariantRecord

if TYPE_CHECKING:
    from .cache import ResultCache

logger = logging.getLogger(__name__)


class VariantScorer:
    """Scores genetic variants using the AlphaGenome API.

    Args:
        config: AlphaVX configuration instance.
    """

    def __init__(self, config: Config) -> None:
        self.config = config
        self._client = None

    def _init_client(self):
        """Lazily initialize the AlphaGenome API client."""
        if self._client is not None:
            return

        from alphagenome.models import dna_client

        self._client = dna_client.create(
            api_key=self.config.api_key,
            address="dns:///gdmscience.googleapis.com:443",
        )
        logger.info("AlphaGenome API client initialized")

    def _get_scorers(self) -> list:
        """Build the list of variant scorers based on configured modalities."""
        from alphagenome.models import variant_scorers

        scorers = []
        for modality in self.config.modalities:
            if modality in variant_scorers.RECOMMENDED_VARIANT_SCORERS:
                scorers.append(variant_scorers.RECOMMENDED_VARIANT_SCORERS[modality])
            else:
                logger.warning("Unknown modality '%s', skipping", modality)
        return scorers

    def score_variant(self, record: VariantRecord) -> pd.DataFrame:
        """Score a single variant across all configured modalities.

        Builds a 1MB genomic interval around the variant, queries the
        AlphaGenome API with all recommended scorers, and returns a
        tidy DataFrame of results.

        Args:
            record: Variant to score.

        Returns:
            DataFrame with columns including biosample_name, gene_name,
            output_type, raw_score, quantile_score, plus chrom/pos/ref/alt.
        """
        from alphagenome.data import genome
        from alphagenome.models import variant_scorers

        self._init_client()

        # Build AlphaGenome objects
        variant = genome.Variant(
            chromosome=record.chrom,
            position=record.pos,
            reference_bases=record.ref,
            alternate_bases=record.alt,
        )
        interval = variant.reference_interval.resize(self.config.sequence_length)

        # Score
        scorers = self._get_scorers()
        logger.info("Scoring %s against %d scorers...", record.key, len(scorers))

        scores_list = self._client.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=scorers,
        )

        # Tidy into DataFrame
        all_dfs = []
        for score_adata in scores_list:
            df = variant_scorers.tidy_scores([score_adata], match_gene_strand=True)
            if df is not None and not df.empty:
                all_dfs.append(df)

        if not all_dfs:
            logger.warning("No scores returned for %s", record.key)
            return pd.DataFrame()

        df = pd.concat(all_dfs, ignore_index=True)

        # Add variant identification columns
        df["chrom"] = record.chrom
        df["pos"] = record.pos
        df["ref"] = record.ref
        df["alt"] = record.alt
        df["variant_key"] = record.key

        return df

    def score_batch(
        self,
        records: list[VariantRecord],
        cache: ResultCache | None = None,
        progress_callback: Callable[[int, int, VariantRecord], None] | None = None,
    ) -> pd.DataFrame:
        """Score a batch of variants with caching and retry logic.

        Iterates through all variants, checks the cache for each,
        and only queries the API for uncached variants. Failed queries
        are retried with exponential backoff. Errors for individual
        variants are logged but don't stop the batch.

        Args:
            records: List of variants to score.
            cache: Optional result cache for skip/store.
            progress_callback: Optional callback(current_index, total, record).

        Returns:
            Combined DataFrame of all scored variants.
        """
        import concurrent.futures
        import threading

        all_dfs: list[pd.DataFrame] = []
        total = len(records)
        cached_count = 0
        error_count = 0

        # Thread-safe counter for progress callback
        completed = 0
        progress_lock = threading.Lock()

        def _process_record(record: VariantRecord) -> pd.DataFrame | None:
            nonlocal cached_count, error_count, completed

            # Check cache first
            if cache is not None and cache.has(record.key):
                cached_data = cache.get(record.key)
                if cached_data is not None:
                    with progress_lock:
                        cached_count += 1
                        completed += 1
                        if progress_callback:
                            progress_callback(completed - 1, total, record)
                    return pd.DataFrame(cached_data)

            # Inform start of scoring if not cached
            with progress_lock:
                if progress_callback:
                    progress_callback(completed, total, record)

            # Score with retry logic
            df = self._score_with_retry(record)

            with progress_lock:
                completed += 1
                if df is not None and not df.empty:
                    # Cache the result
                    if cache is not None:
                        try:
                            cache.put(record.key, df.to_dict(orient="list"))
                        except Exception as e:
                            logger.warning("Failed to cache %s: %s", record.key, e)
                else:
                    error_count += 1

            return df

        # Default max_workers to a sensible limit (e.g. 4-8) to avoid overwhelming the API
        max_workers = min(8, max(4, (os.cpu_count() or 1) + 4))

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Map returns results in the same order as the input
            results = executor.map(_process_record, records)
            for df in results:
                if df is not None and not df.empty:
                    all_dfs.append(df)

        logger.info(
            "Batch complete: %d scored, %d from cache, %d errors out of %d total",
            len(all_dfs) - cached_count,
            cached_count,
            error_count,
            total,
        )

        if not all_dfs:
            return pd.DataFrame()

        return pd.concat(all_dfs, ignore_index=True)

    def _score_with_retry(self, record: VariantRecord) -> pd.DataFrame | None:
        """Score a variant with exponential backoff retries.

        Args:
            record: Variant to score.

        Returns:
            Scoring DataFrame, or None if all retries failed.
        """
        last_error = None
        for attempt in range(self.config.max_retries + 1):
            try:
                return self.score_variant(record)
            except Exception as e:
                last_error = e
                if attempt < self.config.max_retries:
                    delay = self.config.retry_delay * (2**attempt)
                    logger.warning(
                        "Attempt %d/%d failed for %s: %s. Retrying in %.1fs...",
                        attempt + 1,
                        self.config.max_retries + 1,
                        record.key,
                        e,
                        delay,
                    )
                    time.sleep(delay)

        logger.error(
            "All %d attempts failed for %s: %s",
            self.config.max_retries + 1,
            record.key,
            last_error,
        )
        return None
