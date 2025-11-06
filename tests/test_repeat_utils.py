"""Unit tests for repeat utilities."""

import pytest
from strvcf_annotator.core.repeat_utils import (
    extract_repeat_sequence,
    count_repeat_units,
    apply_variant_to_repeat,
    is_perfect_repeat
)


class TestExtractRepeatSequence:
    """Test suite for extract_repeat_sequence."""
    
    def test_simple_repeat(self):
        """Test simple repeat extraction."""
        str_row = {'RU': 'CAG', 'COUNT': 3}
        result = extract_repeat_sequence(str_row)
        assert result == 'CAGCAGCAG'
    
    def test_single_repeat(self):
        """Test single repeat unit."""
        str_row = {'RU': 'AT', 'COUNT': 1}
        result = extract_repeat_sequence(str_row)
        assert result == 'AT'
    
    def test_long_repeat(self):
        """Test longer repeat unit."""
        str_row = {'RU': 'ATCG', 'COUNT': 5}
        result = extract_repeat_sequence(str_row)
        assert result == 'ATCGATCGATCGATCGATCG'


class TestCountRepeatUnits:
    """Test suite for count_repeat_units."""
    
    def test_perfect_repeat(self):
        """Test counting in perfect repeat."""
        assert count_repeat_units('CAGCAGCAG', 'CAG') == 3
    
    def test_imperfect_repeat(self):
        """Test counting in imperfect repeat."""
        assert count_repeat_units('CAGCAGCA', 'CAG') == 2
    
    def test_no_repeats(self):
        """Test sequence with no repeats."""
        assert count_repeat_units('ATCG', 'CAG') == 0
    
    def test_single_repeat(self):
        """Test single repeat unit."""
        assert count_repeat_units('CAG', 'CAG') == 1
    
    def test_overlapping_pattern(self):
        """Test non-overlapping counting."""
        # 'AAA' contains 'AA' once (non-overlapping)
        assert count_repeat_units('AAAA', 'AA') == 2
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        assert count_repeat_units('', 'CAG') == 0
    
    def test_motif_longer_than_sequence(self):
        """Test motif longer than sequence."""
        assert count_repeat_units('CA', 'CAG') == 0


class TestApplyVariantToRepeat:
    """Test suite for apply_variant_to_repeat."""
    
    def test_simple_substitution(self):
        """Test simple substitution."""
        result = apply_variant_to_repeat(100, 'A', 'T', 100, 'AAAA')
        assert result == 'TAAA'
    
    def test_deletion(self):
        """Test deletion."""
        result = apply_variant_to_repeat(100, 'AA', 'A', 100, 'AAAA')
        assert result == 'AAA'
    
    def test_insertion(self):
        """Test insertion."""
        result = apply_variant_to_repeat(100, 'A', 'AT', 100, 'AAAA')
        assert result == 'ATAAA'
    
    def test_variant_in_middle(self):
        """Test variant in middle of repeat."""
        result = apply_variant_to_repeat(102, 'A', 'T', 100, 'AAAA')
        assert result == 'AATA'
    
    def test_variant_before_repeat(self):
        """Test variant starting before repeat."""
        result = apply_variant_to_repeat(98, 'CCA', 'T', 100, 'AAAA')
        assert result == 'AAA'
    
    def test_variant_at_end(self):
        """Test variant at end of repeat."""
        result = apply_variant_to_repeat(103, 'A', 'T', 100, 'AAAA')
        assert result == 'AAAT'


class TestIsPerfectRepeat:
    """Test suite for is_perfect_repeat."""
    
    def test_perfect_repeat(self):
        """Test perfect repeat detection."""
        assert is_perfect_repeat('CAGCAGCAG', 'CAG') is True
    
    def test_imperfect_repeat(self):
        """Test imperfect repeat detection."""
        assert is_perfect_repeat('CAGCAGCA', 'CAG') is False
    
    def test_single_unit(self):
        """Test single repeat unit."""
        assert is_perfect_repeat('CAG', 'CAG') is True
    
    def test_no_repeat(self):
        """Test non-repeat sequence."""
        assert is_perfect_repeat('ATCG', 'CAG') is False
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        assert is_perfect_repeat('', 'CAG') is False
    
    def test_empty_motif(self):
        """Test empty motif."""
        assert is_perfect_repeat('CAG', '') is False
