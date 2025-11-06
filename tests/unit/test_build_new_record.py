"""Unit tests for build_new_record function."""

import pytest
import pysam
import pandas as pd
from strvcf_annotator.core.annotation import build_new_record, make_modified_header
from strvcf_annotator.parsers.generic import GenericParser


@pytest.fixture
def basic_vcf_header():
    """Create a basic VCF header for testing."""
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.contigs.add('chr1', length=1000000)
    header.add_sample('TUMOR')
    header.add_sample('NORMAL')

    # Add basic INFO fields
    header.info.add('DP', 1, 'Integer', 'Total read depth')

    # Add basic FORMAT fields
    header.formats.add('GT', 1, 'String', 'Genotype')
    header.formats.add('AD', 'R', 'Integer', 'Allelic depths')
    header.formats.add('DP', 1, 'Integer', 'Read depth')

    return header


@pytest.fixture
def str_modified_header(basic_vcf_header):
    """Create a mock VCF file and generate STR-modified header."""
    # Create a temporary VCF file to use with make_modified_header
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        # Write header
        f.write(str(basic_vcf_header))
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n')
        # Write a dummy record
        f.write('chr1\t100\t.\tA\tT\t30\tPASS\tDP=50\tGT:AD:DP\t0/1:20,30:50\t0/0:50,0:50\n')
        temp_path = f.name

    try:
        vcf_in = pysam.VariantFile(temp_path)
        header = make_modified_header(vcf_in)
        vcf_in.close()
        return header
    finally:
        os.unlink(temp_path)


@pytest.fixture
def parser():
    """Create a GenericParser instance."""
    return GenericParser()


class TestBuildNewRecordBasic:
    """Basic tests for build_new_record function."""

    def test_returns_variant_record(self, basic_vcf_header, str_modified_header, parser):
        """Test that function returns a VariantRecord."""
        # Create a simple VCF record
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=100,
            alleles=('A', 'T'),
            filter='PASS',
            info={'DP': 50}
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['TUMOR']['AD'] = (20, 30)
        record.samples['TUMOR']['DP'] = 50
        record.samples['NORMAL']['GT'] = (0, 0)
        record.samples['NORMAL']['AD'] = (50, 0)
        record.samples['NORMAL']['DP'] = 50

        # STR metadata
        str_row = {
            'CHROM': 'chr1',
            'START': 95,
            'END': 115,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)
        assert isinstance(new_record, pysam.VariantRecord)

    def test_sets_correct_position(self, basic_vcf_header, str_modified_header, parser):
        """Test that new record has correct position from STR metadata."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=100,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 95,
            'END': 115,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)
        # Position should be START - 1 (0-based)
        assert new_record.start == 94
        assert new_record.stop == 115

    def test_adds_str_info_fields(self, basic_vcf_header, str_modified_header, parser):
        """Test that STR INFO fields are added."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=100,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 95,
            'END': 115,
            'RU': 'CAG',
            'PERIOD': 3,
            'COUNT': 7
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert 'RU' in new_record.info
        assert new_record.info['RU'] == 'CAG'
        assert 'PERIOD' in new_record.info
        assert new_record.info['PERIOD'] == 3
        assert 'REF' in new_record.info
        assert 'PERFECT' in new_record.info

    def test_adds_repcn_format_field(self, basic_vcf_header, str_modified_header, parser):
        """Test that REPCN FORMAT field is added for samples."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=100,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 95,
            'END': 115,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert 'REPCN' in new_record.samples['TUMOR']
        assert 'REPCN' in new_record.samples['NORMAL']
        assert isinstance(new_record.samples['TUMOR']['REPCN'], tuple)


class TestBuildNewRecordAlleles:
    """Tests for allele reconstruction in build_new_record."""

    def test_snv_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test SNV mutation within repeat region."""
        # SNV at position 105
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,  # 0-based
            alleles=('A', 'G')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        # Repeat region: chr1:100-120, repeat unit = "AT"
        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Should have two alleles (reference and mutated repeat sequence)
        assert len(new_record.alleles) == 2
        assert new_record.alleles[0] != new_record.alleles[1]

    def test_insertion_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test insertion mutation."""
        # Insertion at position 105: A -> AAT
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'AAT')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # ALT allele should be longer than REF
        assert len(new_record.alleles[1]) > len(new_record.alleles[0])

    def test_deletion_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test deletion mutation."""
        # Deletion at position 105: AAT -> A
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('AAT', 'A')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # ALT allele should be shorter than REF
        assert len(new_record.alleles[1]) < len(new_record.alleles[0])


class TestBuildNewRecordRepeatCounts:
    """Tests for repeat copy number calculation."""

    def test_calculates_ref_repeat_count(self, basic_vcf_header, str_modified_header, parser):
        """Test that REF INFO field contains repeat count."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert 'REF' in new_record.info
        assert isinstance(new_record.info['REF'], int)
        assert new_record.info['REF'] > 0

    def test_repcn_for_heterozygous(self, basic_vcf_header, str_modified_header, parser):
        """Test REPCN calculation for heterozygous genotype."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'AAT')
        )
        record.samples['TUMOR']['GT'] = (0, 1)  # Het
        record.samples['NORMAL']['GT'] = (0, 0)  # Hom ref

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        tumor_repcn = new_record.samples['TUMOR']['REPCN']
        normal_repcn = new_record.samples['NORMAL']['REPCN']

        # TUMOR should have different copy numbers (0/1 genotype)
        assert len(tumor_repcn) == 2
        assert tumor_repcn[0] != tumor_repcn[1]

        # NORMAL should have same copy numbers (0/0 genotype)
        assert len(normal_repcn) == 2
        assert normal_repcn[0] == normal_repcn[1]

    def test_repcn_for_homozygous_alt(self, basic_vcf_header, str_modified_header, parser):
        """Test REPCN for homozygous ALT genotype."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'AAT')
        )
        record.samples['TUMOR']['GT'] = (1, 1)  # Hom alt
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        tumor_repcn = new_record.samples['TUMOR']['REPCN']

        # Both alleles should have ALT repeat count
        assert tumor_repcn[0] == tumor_repcn[1]


class TestBuildNewRecordPerfectRepeats:
    """Tests for PERFECT field calculation."""

    def test_perfect_field_exists(self, basic_vcf_header, str_modified_header, parser):
        """Test that PERFECT field is always set."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert 'PERFECT' in new_record.info
        assert new_record.info['PERFECT'] in ('TRUE', 'FALSE')


class TestBuildNewRecordEdgeCases:
    """Edge case tests for build_new_record."""

    def test_missing_genotype(self, basic_vcf_header, str_modified_header, parser):
        """Test handling of missing genotype (./.)."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (None, None)
        record.samples['NORMAL']['GT'] = (0, 0)

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Should handle missing genotype gracefully
        assert new_record.samples['TUMOR']['GT'] == (None, None)
        assert new_record.samples['TUMOR']['REPCN'] == (0, 0)

    def test_preserves_original_format_fields(self, basic_vcf_header, str_modified_header, parser):
        """Test that original FORMAT fields are preserved."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['TUMOR']['AD'] = (20, 30)
        record.samples['TUMOR']['DP'] = 50
        record.samples['NORMAL']['GT'] = (0, 0)
        record.samples['NORMAL']['AD'] = (50, 0)
        record.samples['NORMAL']['DP'] = 50

        str_row = {
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        }

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Original FORMAT fields should be preserved
        assert 'AD' in new_record.samples['TUMOR']
        assert 'DP' in new_record.samples['TUMOR']
        assert new_record.samples['TUMOR']['AD'] == (20, 30)
        assert new_record.samples['TUMOR']['DP'] == 50

    def test_with_pandas_series(self, basic_vcf_header, str_modified_header, parser):
        """Test that function works with pandas Series for str_row."""
        record = basic_vcf_header.new_record(
            contig='chr1',
            start=104,
            alleles=('A', 'T')
        )
        record.samples['TUMOR']['GT'] = (0, 1)
        record.samples['NORMAL']['GT'] = (0, 0)

        # Use pandas Series instead of dict
        str_row = pd.Series({
            'CHROM': 'chr1',
            'START': 100,
            'END': 120,
            'RU': 'AT',
            'PERIOD': 2,
            'COUNT': 10
        })

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert isinstance(new_record, pysam.VariantRecord)
        assert new_record.info['RU'] == 'AT'

    def test_different_repeat_units(self, basic_vcf_header, str_modified_header, parser):
        """Test with different repeat unit sizes."""
        test_cases = [
            ('A', 1),      # Mononucleotide
            ('AT', 2),     # Dinucleotide
            ('CAG', 3),    # Trinucleotide
            ('ATCG', 4),   # Tetranucleotide
        ]

        for ru, period in test_cases:
            record = basic_vcf_header.new_record(
                contig='chr1',
                start=104,
                alleles=('A', 'T')
            )
            record.samples['TUMOR']['GT'] = (0, 1)
            record.samples['NORMAL']['GT'] = (0, 0)

            str_row = {
                'CHROM': 'chr1',
                'START': 100,
                'END': 120,
                'RU': ru,
                'PERIOD': period,
                'COUNT': 10
            }

            new_record = build_new_record(record, str_row, str_modified_header, parser)

            assert new_record.info['RU'] == ru
            assert new_record.info['PERIOD'] == period
