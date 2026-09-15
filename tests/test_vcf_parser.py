"""Tests for the VCF parser module."""

import pytest
from pathlib import Path

from alphavx.vcf_parser import parse_vcf, parse_variant_string, VariantRecord


VALID_VCF_CONTENT = """##fileformat=VCFv4.2
##source=test
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr17\t7674220\trs28934578\tG\tA\t100\tPASS\tGENE=TP53
chr13\t32340300\t.\tC\tT\t50\t.\t.
chr7\t117559590\trs113993960\tG\tA,T\t100\tPASS\tGENE=CFTR
"""

MALFORMED_VCF_CONTENT = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT
chr17\tnot_a_number\t.\tG\tA
chr13
"""


class TestParseVcf:
    """Tests for parse_vcf()."""

    def test_parse_valid_vcf(self, tmp_path: Path) -> None:
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(VALID_VCF_CONTENT)

        records = parse_vcf(vcf_file)

        # 3 lines, but chr7 has multi-allelic (G→A and G→T), so 4 records
        assert len(records) == 4

        # First record
        assert records[0].chrom == "chr17"
        assert records[0].pos == 7674220
        assert records[0].ref == "G"
        assert records[0].alt == "A"
        assert records[0].id == "rs28934578"
        assert records[0].info.get("GENE") == "TP53"

    def test_parse_multi_allelic(self, tmp_path: Path) -> None:
        vcf_file = tmp_path / "multi.vcf"
        vcf_file.write_text(VALID_VCF_CONTENT)

        records = parse_vcf(vcf_file)

        # The CFTR line has ALT=A,T → two records at same position
        cftr_records = [r for r in records if r.pos == 117559590]
        assert len(cftr_records) == 2
        assert cftr_records[0].alt == "A"
        assert cftr_records[1].alt == "T"

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            parse_vcf(Path("/nonexistent/file.vcf"))

    def test_empty_vcf_raises(self, tmp_path: Path) -> None:
        vcf_file = tmp_path / "empty.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")

        with pytest.raises(ValueError, match="No valid variant records"):
            parse_vcf(vcf_file)

    def test_malformed_lines_skipped(self, tmp_path: Path) -> None:
        """Malformed lines should be skipped with warnings, not crash."""
        vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT
chr17\tnot_a_number\t.\tG\tA
chr13\t32340300\t.\tC\tT
"""
        vcf_file = tmp_path / "bad.vcf"
        vcf_file.write_text(vcf_content)

        records = parse_vcf(vcf_file)
        # Only the valid line should parse
        assert len(records) == 1
        assert records[0].pos == 32340300

    def test_variant_key_property(self) -> None:
        record = VariantRecord(chrom="chr1", pos=100, ref="A", alt="G")
        assert record.key == "chr1:100:A>G"

    def test_parse_gzipped_vcf(self, tmp_path: Path) -> None:
        """Should correctly parse .vcf.gz files."""
        import gzip

        vcf_file = tmp_path / "test.vcf.gz"
        with gzip.open(vcf_file, "wt") as f:
            f.write(VALID_VCF_CONTENT)

        records = parse_vcf(vcf_file)
        assert len(records) == 4
        assert records[0].chrom == "chr17"
        assert records[0].pos == 7674220


class TestParseVariantString:
    """Tests for parse_variant_string()."""

    def test_arrow_format(self) -> None:
        record = parse_variant_string("chr17:7674220:G>A")
        assert record.chrom == "chr17"
        assert record.pos == 7674220
        assert record.ref == "G"
        assert record.alt == "A"

    def test_colon_format(self) -> None:
        record = parse_variant_string("chr17:7674220:G:A")
        assert record.chrom == "chr17"
        assert record.pos == 7674220
        assert record.ref == "G"
        assert record.alt == "A"

    def test_whitespace_stripped(self) -> None:
        record = parse_variant_string("  chr1:500:C>T  ")
        assert record.chrom == "chr1"
        assert record.pos == 500

    def test_invalid_format_raises(self) -> None:
        with pytest.raises(ValueError, match="Cannot parse variant"):
            parse_variant_string("this-is-not-a-variant")

    def test_invalid_position_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid position"):
            parse_variant_string("chr1:abc:G>A")

    def test_missing_alt_raises(self) -> None:
        with pytest.raises(ValueError):
            parse_variant_string("chr1:100:G>")
