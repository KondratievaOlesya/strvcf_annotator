# API Documentation

## Overview

The `strvcf_annotator` package provides both a library API and command-line interface for annotating VCF files with STR (Short Tandem Repeat) information.

## Library API

### STRAnnotator Class

The main class for STR annotation functionality.

```python
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator(str_bed_path, parser=None)
```

#### Parameters

- `str_bed_path` (str): Path to BED file containing STR regions
- `parser` (BaseVCFParser, optional): Custom parser for genotype extraction. Uses GenericParser if None.

#### Methods

##### annotate_vcf_file(input_path, output_path)

Annotate a single VCF file.

```python
annotator.annotate_vcf_file('input.vcf', 'output.vcf')
```

**Parameters:**
- `input_path` (str): Path to input VCF file
- `output_path` (str): Path to output VCF file

**Example:**
```python
annotator = STRAnnotator('repeats.bed')
annotator.annotate_vcf_file('input.vcf', 'output.vcf')
```

##### annotate_vcf_stream(vcf_in)

Generator that yields annotated VCF records from a stream.

```python
import pysam

vcf_in = pysam.VariantFile('input.vcf')
for record in annotator.annotate_vcf_stream(vcf_in):
    # Process record
    pass
```

**Parameters:**
- `vcf_in` (pysam.VariantFile): VCF file object

**Yields:**
- `pysam.VariantRecord`: Annotated VCF records

**Example:**
```python
import pysam

annotator = STRAnnotator('repeats.bed')
vcf_in = pysam.VariantFile('input.vcf')

for record in annotator.annotate_vcf_stream(vcf_in):
    print(f"Repeat unit: {record.info['RU']}")
    print(f"Copy number: {record.info['REF']}")
```

##### process_directory(input_dir, output_dir)

Batch process all VCF files in a directory.

```python
annotator.process_directory('vcf_files/', 'annotated_vcfs/')
```

**Parameters:**
- `input_dir` (str): Directory containing input VCF files
- `output_dir` (str): Directory for output VCF files

**Example:**
```python
annotator = STRAnnotator('repeats.bed')
annotator.process_directory('vcf_files/', 'annotated/')
```

##### get_str_at_position(chrom, pos)

Get STR region at a specific genomic position.

```python
str_region = annotator.get_str_at_position(chrom, pos)
```

**Parameters:**
- `chrom` (str): Chromosome name
- `pos` (int): Genomic position (1-based)

**Returns:**
- `dict` or `None`: STR region data if position is within an STR

**Example:**
```python
annotator = STRAnnotator('repeats.bed')
str_region = annotator.get_str_at_position('chr1', 1000000)

if str_region:
    print(f"Repeat unit: {str_region['RU']}")
    print(f"Period: {str_region['PERIOD']}")
```

##### get_statistics()

Get statistics about loaded STR regions.

```python
stats = annotator.get_statistics()
```

**Returns:**
- `dict`: Statistics including total regions, chromosomes, repeat units

**Example:**
```python
annotator = STRAnnotator('repeats.bed')
stats = annotator.get_statistics()

print(f"Всего регионов: {stats['total_regions']}")
print(f"Хромосом: {stats['chromosomes']}")
print(f"Уникальных единиц повтора: {stats['unique_repeat_units']}")
print(f"Среднее количество повторов: {stats['mean_repeat_count']:.2f}")
```


### Utility Function

#### annotate_vcf(input_vcf, str_bed, output_vcf, parser=None)

Simple function interface for annotating a single VCF file.

```python
from strvcf_annotator import annotate_vcf

annotate_vcf('input.vcf', 'str_regions.bed', 'output.vcf', parser=None)
```

**Parameters:**
- `input_vcf` (str): Path to input VCF file
- `str_bed` (str): Path to STR BED file
- `output_vcf` (str): Path to output VCF file
- `parser` (BaseVCFParser, optional): Custom parser

## Parser Interface

### BaseVCFParser

Abstract base class for VCF parsers.

```python
from strvcf_annotator.parsers.base import BaseVCFParser

class CustomParser(BaseVCFParser):
    def get_genotype(self, record, sample_idx):
        # Implementation
        pass
    
    def has_variant(self, record, sample_idx):
        # Implementation
        pass
    
    def extract_info(self, record, sample_idx):
        # Implementation
        pass
    
    def validate_record(self, record):
        # Implementation
        pass
```

#### Methods to Implement

##### get_genotype()

Extracts genotype from a VCF record.

```python
def get_genotype(self, record, sample_idx):
    """
    Parameters:
        record: pysam.VariantRecord
        sample_idx: int - sample index
    
    Returns:
        Optional[Tuple[int, int]] - genotype as (allele1, allele2) or None
    """
```

##### has_variant()

Checks if the variant is present in the sample.

```python
def has_variant(self, record, sample_idx):
    """
    Parameters:
        record: pysam.VariantRecord
        sample_idx: int - sample index
    
    Returns:
        bool - True if the variant is present
    """
```

##### extract_info()

Extracts additional FORMAT fields.

```python
def extract_info(self, record, sample_idx):
    """
    Parameters:
        record: pysam.VariantRecord
        sample_idx: int - sample index
    
    Returns:
        Dict[str, Any] - dictionary with FORMAT fields
    """
```

##### validate_record()

Checks if the record is compatible with the parser.

```python
def validate_record(self, record):
    """
    Parameters:
        record: pysam.VariantRecord
    
    Returns:
        bool - True if the record is valid
    """
```

### GenericParser

Built-in parser for standard VCF format fields (GT, AD, DP).

```python
from strvcf_annotator.parsers.generic import GenericParser

parser = GenericParser()
annotator = STRAnnotator('str_regions.bed', parser=parser)
```

## Core Functions

### STR Reference

#### load_str_reference(str_path)

Load STR reference data from BED file.

```python
from strvcf_annotator.core.str_reference import load_str_reference

str_df = load_str_reference('str_regions.bed')
```

**Parameters:**
- `str_path` (str): Path to BED file

**Returns:**
- `pd.DataFrame`: DataFrame with STR regions

### Repeat Utilities

#### extract_repeat_sequence(str_row)

Reconstruct repeat sequence from STR metadata.

```python
from strvcf_annotator.core.repeat_utils import extract_repeat_sequence

str_row = {'RU': 'CAG', 'COUNT': 5}
seq = extract_repeat_sequence(str_row)  # Returns 'CAGCAGCAGCAGCAG'
```

#### count_repeat_units(sequence, motif)

Count non-overlapping occurrences of repeat motif.

```python
from strvcf_annotator.core.repeat_utils import count_repeat_units

count = count_repeat_units('CAGCAGCAG', 'CAG')  # Returns 3
```

#### apply_variant_to_repeat(pos, ref, alt, repeat_start, repeat_seq)

Apply variant to repeat sequence.

```python
from strvcf_annotator.core.repeat_utils import apply_variant_to_repeat

mutated = apply_variant_to_repeat(100, 'A', 'T', 100, 'AAAA')  # Returns 'TAAA'
```

#### is_perfect_repeat(sequence, motif)

Check if sequence is a perfect repeat of the motif.

```python
from strvcf_annotator.core.repeat_utils import is_perfect_repeat

is_perfect = is_perfect_repeat('CAGCAGCAG', 'CAG')  # Returns True
```

## Exceptions

### ValidationError

Raised for validation errors.

```python
from strvcf_annotator.utils.validation import ValidationError

try:
    annotator = STRAnnotator('nonexistent.bed')
except ValidationError as e:
    print(f"Validation error: {e}")
```

### Validation Functions

```python
from strvcf_annotator.utils.validation import (
    validate_file_path,
    validate_vcf_file,
    validate_str_bed_file
)

# Validate file path
path = validate_file_path('input.vcf', must_exist=True)

# Validate VCF file
is_valid = validate_vcf_file('input.vcf')

# Validate STR BED file
is_valid = validate_str_bed_file('repeats.bed')
```


## Output Format

The annotated VCF includes the following additional fields:

### INFO Fields

- `RU` (String): Repeat unit sequence
- `PERIOD` (Integer): Repeat period (length of unit)
- `REF` (Integer): Reference copy number
- `PERFECT` (String): "TRUE" if both REF and ALT are perfect repeats, "FALSE" otherwise

### FORMAT Fields

- `REPCN` (Integer, 2 values): Genotype as number of repeat motif copies
