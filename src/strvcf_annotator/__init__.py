"""Top-level package for strvcf_annotator.

STR (Short Tandem Repeat) annotation tool for VCF files.
Provides both library and CLI interfaces for annotating variants
that overlap with STR regions.
"""

__author__ = """Olesia Kondrateva"""
__email__ = 'xkdnoa@gmail.com'
__version__ = '0.1.0'

# Public API exports
from .api import STRAnnotator, annotate_vcf
from .parsers.base import BaseVCFParser
from .parsers.generic import GenericParser
from .core.str_reference import load_str_reference
from .core.repeat_utils import (
    extract_repeat_sequence,
    count_repeat_units,
    apply_variant_to_repeat,
    is_perfect_repeat
)
from .utils.validation import ValidationError

__all__ = [
    # Main API
    'STRAnnotator',
    'annotate_vcf',
    
    # Parsers
    'BaseVCFParser',
    'GenericParser',
    
    # Core functions
    'load_str_reference',
    'extract_repeat_sequence',
    'count_repeat_units',
    'apply_variant_to_repeat',
    'is_perfect_repeat',
    
    # Exceptions
    'ValidationError',
]
