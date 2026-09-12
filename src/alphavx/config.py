"""Configuration loading and validation for AlphaVX."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

# Load .env.alphavx file if present (before any os.environ reads)
# NOTE: We use ".env.alphavx" instead of ".env" because anndata (an
# alphagenome dependency) uses pydantic-settings which auto-reads ".env"
# and rejects unknown keys like ALPHAVX_API_KEY.
try:
    from dotenv import load_dotenv

    _ENV_FILENAME = ".env.alphavx"

    # Search current dir and project root
    _candidates = [
        Path.cwd() / _ENV_FILENAME,
        Path(__file__).resolve().parent.parent.parent / _ENV_FILENAME,
    ]
    for _candidate in _candidates:
        if _candidate.exists():
            load_dotenv(_candidate)
            break
except ImportError:
    pass

logger = logging.getLogger(__name__)

DEFAULT_MODALITIES = [
    "RNA_SEQ",
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
    "DNASE",
    "ATAC",
    "CHIP_HISTONE",
    "CHIP_TF",
]


@dataclass
class Config:
    """AlphaVX runtime configuration."""

    api_key: str = ""
    max_retries: int = 3
    retry_delay: float = 5.0
    sequence_length: int = 2**20  # 1,048,576 — AlphaGenome optimal
    quantile_threshold: float = 0.995
    output_dir: Path = field(default_factory=lambda: Path("results"))
    cache_enabled: bool = True
    modalities: list[str] = field(default_factory=lambda: list(DEFAULT_MODALITIES))


def load_config(config_path: Path | None = None) -> Config:
    """Load configuration from YAML file and/or environment variables.

    Priority: YAML file values > environment variables > defaults.

    Args:
        config_path: Optional path to an alphavx.yaml configuration file.

    Returns:
        Populated Config instance.

    Raises:
        ValueError: If no API key is found in config or environment.
    """
    config = Config()

    # Load from YAML if provided
    if config_path is not None:
        try:
            import yaml
        except ImportError:
            logger.warning("pyyaml not installed — ignoring config file %s", config_path)
        else:
            with open(config_path) as f:
                raw = yaml.safe_load(f) or {}

            api_section = raw.get("api", {})
            scoring_section = raw.get("scoring", {})
            output_section = raw.get("output", {})

            if "max_retries" in api_section:
                config.max_retries = int(api_section["max_retries"])
            if "retry_delay" in api_section:
                config.retry_delay = float(api_section["retry_delay"])
            if "modalities" in scoring_section:
                config.modalities = scoring_section["modalities"]
            if "quantile_threshold" in scoring_section:
                config.quantile_threshold = float(scoring_section["quantile_threshold"])
            if "sequence_length" in scoring_section:
                config.sequence_length = int(scoring_section["sequence_length"])
            if "cache" in output_section:
                config.cache_enabled = bool(output_section["cache"])

            logger.info("Loaded config from %s", config_path)

    # Resolve API key: env var takes precedence if config doesn't set it
    if not config.api_key:
        config.api_key = os.environ.get(
            "ALPHAVX_API_KEY",
            os.environ.get("ALPHAGENOME_API_KEY", ""),
        )

    if not config.api_key:
        raise ValueError(
            "No API key found. Set ALPHAVX_API_KEY or ALPHAGENOME_API_KEY environment "
            "variable, or provide it in your alphavx.yaml config file.\n"
            "Sign up at: https://deepmind.google.com/science/alphagenome/"
        )

    return config
