"""Core modules for STR annotation functionality."""

from .str_reference import load_str_reference
from .repeat_utils import extract_repeat_sequence, count_repeat_units, apply_variant_to_repeat
from .annotation import make_modified_header, build_new_record, should_skip_genotype
from .vcf_processor import (
    check_vcf_sorted,
    reset_and_sort_vcf,
    generate_annotated_records,
    annotate_vcf_to_file
)

__all__ = [
    'load_str_reference',
    'extract_repeat_sequence',
    'count_repeat_units',
    'apply_variant_to_repeat',
    'make_modified_header',
    'build_new_record',
    'should_skip_genotype',
    'check_vcf_sorted',
    'reset_and_sort_vcf',
    'generate_annotated_records',
    'annotate_vcf_to_file'
]
