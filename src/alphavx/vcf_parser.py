"""VCF file parsing and variant string parsing for AlphaVX."""

from __future__ import annotations

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class VariantRecord:
    """A single genetic variant parsed from VCF or command-line input."""

    chrom: str
    pos: int
    ref: str
    alt: str
    id: str | None = None
    quality: float | None = None
    filter: str | None = None
    info: dict[str, str] = field(default_factory=dict)

    @property
    def key(self) -> str:
        """Unique string key for this variant: chr:pos:ref>alt."""
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"

    def __str__(self) -> str:
        return self.key


def parse_vcf(vcf_path: Path) -> list[VariantRecord]:
    """Parse a VCF file into a list of VariantRecord objects.

    Handles multi-allelic sites by splitting into separate records.
    Skips malformed lines with a warning rather than crashing.

    Args:
        vcf_path: Path to a VCF file (.vcf or .vcf.gz).

    Returns:
        List of parsed VariantRecord objects.

    Raises:
        FileNotFoundError: If vcf_path does not exist.
        ValueError: If the file contains no valid variant lines.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    records: list[VariantRecord] = []

    opener = gzip.open(vcf_path, "rt") if str(vcf_path).endswith(".gz") else open(vcf_path)
    with opener as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()

            # Skip metadata and header lines
            if line.startswith("##") or line.startswith("#CHROM") or not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                logger.warning("Line %d: too few columns (%d), skipping", line_num, len(parts))
                continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                logger.warning("Line %d: invalid POS '%s', skipping", line_num, parts[1])
                continue

            variant_id = parts[2] if parts[2] != "." else None
            ref = parts[3]
            alt_field = parts[4]

            # Parse QUAL
            quality = None
            if len(parts) > 5 and parts[5] != ".":
                try:
                    quality = float(parts[5])
                except ValueError:
                    pass

            # Parse FILTER
            filt = parts[6] if len(parts) > 6 and parts[6] != "." else None

            # Parse INFO into a dict
            info: dict[str, str] = {}
            if len(parts) > 7 and parts[7] != ".":
                for item in parts[7].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info[k] = v
                    else:
                        info[item] = ""

            # Handle multi-allelic: split comma-separated ALTs into separate records
            for alt in alt_field.split(","):
                alt = alt.strip()
                if not alt or alt == ".":
                    continue
                records.append(
                    VariantRecord(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        id=variant_id,
                        quality=quality,
                        filter=filt,
                        info=dict(info),
                    )
                )

    if not records:
        raise ValueError(f"No valid variant records found in {vcf_path}")

    logger.info("Parsed %d variants from %s", len(records), vcf_path)
    return records


def parse_variant_string(variant_str: str) -> VariantRecord:
    """Parse a variant string like 'chr17:7674220:G>A' or 'chr17:7674220:G:A'.

    Accepted formats:
        - chr:pos:ref>alt  (e.g., chr17:7674220:G>A)
        - chr:pos:ref:alt  (e.g., chr17:7674220:G:A)

    Args:
        variant_str: Variant specification string.

    Returns:
        Parsed VariantRecord.

    Raises:
        ValueError: If the string cannot be parsed.
    """
    variant_str = variant_str.strip()

    # Try chr:pos:ref>alt format
    if ">" in variant_str:
        parts = variant_str.split(":")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid variant format: '{variant_str}'. "
                "Expected chr:pos:ref>alt (e.g., chr17:7674220:G>A)"
            )
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            raise ValueError(f"Invalid position: '{parts[1]}' in variant '{variant_str}'")

        ref_alt = parts[2].split(">")
        if len(ref_alt) != 2 or not ref_alt[0] or not ref_alt[1]:
            raise ValueError(
                f"Invalid ref>alt: '{parts[2]}' in variant '{variant_str}'. "
                "Expected format like G>A"
            )
        return VariantRecord(chrom=chrom, pos=pos, ref=ref_alt[0], alt=ref_alt[1])

    # Try chr:pos:ref:alt format
    parts = variant_str.split(":")
    if len(parts) == 4:
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            raise ValueError(f"Invalid position: '{parts[1]}' in variant '{variant_str}'")
        ref, alt = parts[2], parts[3]
        if not ref or not alt:
            raise ValueError(f"Empty ref or alt in variant '{variant_str}'")
        return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt)

    raise ValueError(
        f"Cannot parse variant: '{variant_str}'. "
        "Expected chr:pos:ref>alt or chr:pos:ref:alt "
        "(e.g., chr17:7674220:G>A or chr17:7674220:G:A)"
    )
