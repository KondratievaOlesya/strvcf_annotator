"""Utility functions for VCF processing."""

import pysam
from typing import Dict, Any


def normalize_info_fields(record: pysam.VariantRecord, header: pysam.VariantHeader) -> Dict[str, Any]:
    """Normalize INFO fields for proper VCF serialization.
    
    Handles various INFO field types and ensures they are properly formatted
    for writing to VCF files. Handles Flags, Strings, and R-type fields specially.
    
    Parameters
    ----------
    record : pysam.VariantRecord
        VCF record with INFO fields to normalize
    header : pysam.VariantHeader
        VCF header with field definitions
        
    Returns
    -------
    Dict[str, Any]
        Normalized INFO fields ready for VCF writing
        
    Notes
    -----
    - Flag fields are included only if True
    - String fields with multiple values are joined with "|"
    - R-type fields (REF + ALT) are clipped to first 2 values
    - Unknown fields are skipped
    """
    fixed_info = {}

    for key, val in record.info.items():
        if key not in header.info:
            continue  # Skip unknown fields

        desc = header.info[key]

        # Handle Flags (should be included as key only if present)
        if desc.type == "Flag":
            if val:
                fixed_info[key] = True
            continue

        # Handle single String field (but given as tuple or list)
        if desc.type == "String" and desc.number == 1:
            if isinstance(val, (list, tuple)):
                fixed_info[key] = "|".join(map(str, val))
            else:
                fixed_info[key] = val
            continue
        
        # Handle single Integer/Float field (but given as tuple or list)
        if desc.number == 1 and isinstance(val, (list, tuple)):
            fixed_info[key] = val[0] if len(val) > 0 else 0
            continue
        
        # Clip R fields (REF + ALT)
        if desc.number == "R" and isinstance(val, (list, tuple)):
            fixed_info[key] = list(val[:2])
            continue

        # Pass through other values
        fixed_info[key] = val

    return fixed_info


def get_sample_by_name(record: pysam.VariantRecord, sample_name: str) -> Any:
    """Get sample data by name from VCF record.
    
    Parameters
    ----------
    record : pysam.VariantRecord
        VCF record
    sample_name : str
        Name of sample to retrieve
        
    Returns
    -------
    Any
        Sample data object
        
    Raises
    ------
    KeyError
        If sample name not found in record
    """
    return record.samples[sample_name]


def get_sample_by_index(record: pysam.VariantRecord, sample_idx: int) -> Any:
    """Get sample data by index from VCF record.
    
    Parameters
    ----------
    record : pysam.VariantRecord
        VCF record
    sample_idx : int
        Index of sample to retrieve
        
    Returns
    -------
    Any
        Sample data object
        
    Raises
    ------
    IndexError
        If sample index out of range
    """
    samples = list(record.samples.values())
    return samples[sample_idx]


def has_format_field(record: pysam.VariantRecord, field_name: str) -> bool:
    """Check if FORMAT field exists in any sample.
    
    Parameters
    ----------
    record : pysam.VariantRecord
        VCF record
    field_name : str
        Name of FORMAT field to check
        
    Returns
    -------
    bool
        True if field exists in at least one sample
    """
    for sample in record.samples.values():
        if field_name in sample:
            return True
    return False
