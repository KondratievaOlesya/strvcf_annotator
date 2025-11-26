"""Unit tests for VCF parsers."""

import pysam
import pytest

from strvcf_annotator.parsers.generic import GenericParser


class TestGenericParser:
    """Test suite for GenericParser."""

    @pytest.fixture
    def parser(self):
        """Create GenericParser instance."""
        return GenericParser()

    @pytest.fixture
    def vcf_header(self):
        """Create minimal VCF header."""
        header = pysam.VariantHeader()
        header.add_sample("Sample1")
        header.add_sample("Sample2")
        header.formats.add("GT", 1, "String", "Genotype")
        header.formats.add("AD", "R", "Integer", "Allelic depth")
        header.formats.add("DP", 1, "Integer", "Total depth")
        header.contigs.add("chr1")
        return header

    def test_get_genotype_valid(self, parser, vcf_header):
        """Test genotype extraction with valid GT."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (0, 1)

        gt = parser.get_genotype(record, 0)
        assert gt == (0, 1)

    def test_get_genotype_missing(self, parser, vcf_header):
        """Test genotype extraction with missing GT."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (None, None)

        gt = parser.get_genotype(record, 0)
        assert gt is None

    def test_get_genotype_invalid_index(self, parser, vcf_header):
        """Test genotype extraction with invalid sample index."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (0, 1)

        gt = parser.get_genotype(record, 10)
        assert gt is None

    def test_has_variant_with_alt(self, parser, vcf_header):
        """Test variant detection with ALT allele."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (0, 1)

        assert parser.has_variant(record, 0) is True

    def test_has_variant_homozygous_ref(self, parser, vcf_header):
        """Test variant detection with homozygous reference."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (0, 0)

        assert parser.has_variant(record, 0) is False

    def test_has_variant_with_ad(self, parser, vcf_header):
        """Test variant detection using AD when GT missing."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["GT"] = (None, None)
        record.samples["Sample1"]["AD"] = (10, 5)

        assert parser.has_variant(record, 0) is True

    def test_extract_info(self, parser, vcf_header):
        """Test extraction of FORMAT fields."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["Sample1"]["AD"] = (10, 5)
        record.samples["Sample1"]["DP"] = 15

        info = parser.extract_info(record, 0)
        assert info["AD"] == (10, 5)
        assert info["DP"] == 15

    def test_validate_record_valid(self, parser, vcf_header):
        """Test record validation with valid record."""
        record = vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))

        assert parser.validate_record(record) is True

    def test_validate_record_no_samples(self, parser):
        """Test record validation with no samples."""
        header = pysam.VariantHeader()
        header.contigs.add("chr1")
        record = header.new_record(contig="chr1", start=100, alleles=("A", "T"))

        assert parser.validate_record(record) is False

    def test_validate_record_none(self, parser):
        """Test record validation with None."""
        assert parser.validate_record(None) is False
