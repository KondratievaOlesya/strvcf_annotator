# Testing Documentation

## Overview

The `strvcf_annotator` project has a comprehensive test suite with 158 tests across multiple categories, achieving 98.7% test coverage (156 passing, 2 known pre-existing issues).

**Current Test Statistics:**
- **Total Tests**: 158
- **Passing**: 156 (98.7%)
- **Known Issues**: 2 (pre-existing coordinate system bugs)
- **Test Files**: 15 files 
- **Test Data**: 63MB across 7 VCF/BED files
- **Categories**: Unit, Integration, CLI, Performance

## Quick Start

```bash
# Run all tests
pytest tests/

# Run with coverage report
pytest tests/ --cov=src/strvcf_annotator --cov-report=html

# Run specific test category
pytest tests/unit/              # Unit tests only
pytest tests/integration/       # Integration tests only
pytest tests/cli/              # CLI tests only

# Run specific test file
pytest tests/unit/test_annotation.py

# Run specific test
pytest tests/unit/test_annotation.py::TestBuildNewRecord::test_creates_new_record

# Run with verbose output
pytest tests/ -v

# Run and stop at first failure
pytest tests/ -x

# Update VCF hash values (after intentional changes)
pytest tests/integration/test_vcf_pipeline.py --update-vcf-hashes
```

## Test Structure

```
tests/
├── conftest.py                          # Global pytest configuration and fixtures
├── data/                                # Test data files (62.1MB)
│   ├── GRCh38_repeats.bed               # STR reference panel (10,615 repeats)
│   ├── test.vcf.gz                      # TCGA-COAD VCF (7,672 variants)
│   ├── pindel_header.vcf                # Pindel VCF header (1 variant)
│   ├── mutec2_indel.vcf.gz              # Mutect2 indels with index (40,016 variants)
│   └── TCGA-DC-6682.vcf                 # TCGA annotated with GangSTR VCF (181,563 variants)
├── output/                              # Generated test output files
├── unit/                                # Unit tests (isolated component testing)
│   ├── test_annotation.py              # 21 tests - build_new_record, make_modified_header, should_skip_genotype
│   ├── test_build_new_record.py        # 17 tests - record building with various edge cases
│   ├── test_make_modified_header.py    # 6 tests - VCF header modification
│   ├── test_repeat_utils.py            # 34 tests - repeat sequence manipulation
│   ├── test_str_reference.py           # 8 tests - BED file loading and validation
│   ├── test_validation.py              # 13 tests - input file validation
│   └── test_vcf_processor.py           # 27 tests - VCF processing functions
├── integration/                         # End-to-end tests
│   └── test_vcf_pipeline.py            # 4 tests - full pipeline with hash validation
├── cli/                                 # Command-line interface tests
│   └── test_cli.py                     # 15 tests - CLI argument parsing and execution
├── parsers/                             # Parser-specific tests
│   └── test_parsers.py                 # 6 tests - GenericParser and MutectParser
├── performance/                         # Performance benchmarks
│   └── test_performance.py             # 5 tests - timing benchmarks for large files
└── test_integration.py                  # 2 tests - API usage patterns

Total: 158 tests across 15 files
```

## Test Data Files

| File | Type | Variants | Caller | Purpose |
|------|------|----------|--------|---------|
| `GRCh38_repeats.bed` | BED | 10,615 STRs | - | Reference panel for GRCh38 |
| `test.vcf.gz` | VCF.gz | 7,672 | Unknown | Somatic SNPs (compressed) |
| `pindel_header.vcf` | VCF | 1 | Pindel | Pindel header for future parsers |
| `mutec2_indel.vcf.gz` | VCF.gz | 40,016 | Mutect2 | With .tbi index file |
| `TCGA-DC-6682.vcf` | VCF | 181,563 | GangSTR | TCGA annotated with GangSTR VCF |

**Total Test Data Size**: ~62 MB

## Running Tests

### By Category

```bash
# Unit tests - Fast, isolated component testing
pytest tests/unit/ -v

# Integration tests - End-to-end pipeline validation
pytest tests/integration/ -v

# CLI tests - Command-line interface testing
pytest tests/cli/ -v

# Parser tests - VCF parser validation
pytest tests/parsers/ -v

# Performance tests - Benchmarking (may take several minutes)
pytest tests/performance/ -v
```

### By Module

```bash
# Annotation engine tests
pytest tests/unit/test_annotation.py

# Repeat manipulation tests
pytest tests/unit/test_repeat_utils.py

# VCF processor tests
pytest tests/unit/test_vcf_processor.py

# STR reference loading tests
pytest tests/unit/test_str_reference.py

# Input validation tests
pytest tests/unit/test_validation.py

# Full pipeline tests (with hash validation)
pytest tests/integration/test_vcf_pipeline.py
```

### With Coverage

```bash
# Generate HTML coverage report
pytest tests/ --cov=src/strvcf_annotator --cov-report=html
# Open htmlcov/index.html to view detailed coverage

# Terminal coverage report
pytest tests/ --cov=src/strvcf_annotator --cov-report=term-missing

# Coverage for specific module
pytest tests/unit/test_annotation.py --cov=src/strvcf_annotator.core.annotation
```

### Debugging Failed Tests

```bash
# Show detailed output including print statements
pytest tests/ -v -s

# Stop at first failure
pytest tests/ -x

# Show local variables on failure
pytest tests/ -l

# Run only failed tests from last run
pytest tests/ --lf

# Debug with pdb on failure
pytest tests/ --pdb
```

## Test Scenarios

### CLI Testing (`test_cli.py`)

Tests command-line argument parsing and execution:

- **Basic execution**: `strvcf-annotator --input in.vcf --str-bed repeats.bed --output out.vcf`
- **Logging**: `--log-level DEBUG/INFO/WARNING/ERROR`
- **File validation**: Missing files, invalid formats
- **Help and version**: `--help`, `--version`

### VCF Processing (`test_vcf_processor.py`)

Tests core VCF processing functions:

- **Sorting detection**: `check_vcf_sorted()` - identifies unsorted VCFs
- **VCF sorting**: `reset_and_sort_vcf()` - sorts records by contig order
- **Record generation**: `generate_annotated_records()` - filters and annotates variants
- **Edge cases**: Empty VCFs, single records, multi-chromosome files

### Variant Annotation (`test_annotation.py`, `test_build_new_record.py`)

Tests STR annotation logic:

- **Header modification**: `make_modified_header()` - adds RU, PERIOD, REF, PERFECT INFO fields
- **Record building**: `build_new_record()` - creates annotated records with full repeat sequences
- **Genotype filtering**: `should_skip_genotype()` - somatic mode filtering logic
- **Edge cases**:
  - Deletions (GAA>G removes 2bp)
  - Insertions (A>AAT adds 2bp)
  - Substitutions (G>C point mutations)
  - Multi-sample VCFs
  - Missing genotypes

### Genotype Handling

Tests genotype extraction and REPCN calculation:

- **Generic parser**: Handles standard VCF GT format
- **Mutect parser**: Handles Mutect2-specific FORMAT fields
- **Somatic filtering**: Skips variants where tumor==normal genotypes
- **Heterozygous variants**: Correctly calculates (REF_count, ALT_count)
- **Homozygous variants**: Correctly calculates (ALT_count, ALT_count)

### Repeat Sequence Manipulation (`test_repeat_utils.py`)

Tests repeat sequence utilities:

- **Sequence extraction**: `extract_repeat_sequence()` - gets repeat from BED coordinates
- **Variant application**: `apply_variant_to_repeat()` - applies VCF mutation to repeat
- **Repeat counting**: `count_repeat_units()` - counts RU copies (e.g., "ATATATAT" = 4× "AT")
- **Perfect repeat detection**: `is_perfect_repeat()` - validates exact RU repetition
- **Edge cases**: Partial repeats, complex insertions/deletions

### STR Reference Loading (`test_str_reference.py`)

Tests BED file loading and validation:

- **File loading**: `load_str_reference()` - parses BED into DataFrame
- **Column validation**: Ensures CHROM, START, END, PERIOD, RU, COUNT columns
- **Coordinate validation**: START < END, positive coordinates
- **Reference sequence validation**: Ensures reference FASTA access

## Pytest Fixtures

### Global Fixtures (`conftest.py`)

These fixtures are available to all tests:

```python
@pytest.fixture(scope='session')
def data_dir():
    """Path to tests/data/ directory."""

@pytest.fixture(scope='session')
def output_dir():
    """Path to tests/output/ directory (created if needed)."""

@pytest.fixture
def sample_str_bed(tmp_path):
    """Creates temporary BED file with 3 STR regions."""

@pytest.fixture
def sample_vcf(tmp_path):
    """Creates temporary VCF file with 3 variants."""
```

### Common Unit Test Fixtures

```python
@pytest.fixture
def basic_vcf_header():
    """Creates basic VCF header with chr1/chr2/chr3, TUMOR/NORMAL samples."""

@pytest.fixture
def sorted_vcf_file(basic_vcf_header):
    """Creates temporary sorted VCF file using pysam API."""

@pytest.fixture
def unsorted_vcf_file(basic_vcf_header):
    """Creates temporary unsorted VCF file."""

@pytest.fixture
def str_dataframe():
    """Creates sample STR DataFrame with 3 regions."""
```

### Parametrized Fixtures

```python
@pytest.fixture(params=[...])
def vcf_case(request, data_dir):
    """Parametrized fixture - runs test 4 times with different VCF files."""
    # Used in test_vcf_pipeline.py for hash validation

@pytest.fixture(params=[...])
def vcf_in(request, data_dir):
    """Parametrized fixture - runs test with 4 different VCF formats."""
    # Used in test_make_modified_header.py
```

## Performance Benchmarks

Performance tests validate that the tool can process large files efficiently:

```python
# test_performance.py benchmarks

def test_annotation_performance_large_vcf():
    """Benchmark: Process TCGA VCF (181,563 variants) in < 60 seconds."""

def test_annotation_performance_compressed():
    """Benchmark: Process compressed VCF.gz efficiently."""

def test_memory_usage():
    """Benchmark: Memory usage stays reasonable for large files."""
```

**Baseline Performance** (on reference hardware):
- Large VCF (181K variants): ~30-45 seconds
- Compressed VCF.gz: ~10-15 seconds
- Memory usage: <500 MB for large files

### Guidelines

1. **Use descriptive test names**: `test_filters_variants_with_identical_genotypes()` not `test_filter()`
2. **One assertion per concept**: Test one thing at a time for clarity
3. **Use fixtures for setup**: Avoid duplicating setup code across tests
4. **Test edge cases**: Empty inputs, single items, large datasets
5. **Test error conditions**: Invalid inputs should raise appropriate errors
6. **Use pysam API for VCF files**: Don't manipulate VCF text directly
7. **Clean up resources**: Use `tmp_path` fixture for temporary files
8. **Document test purpose**: Add docstrings explaining what is tested

## Hash-Based Validation

Integration tests use MD5 hashes to validate output correctness:

### Current Hash Values (`test_vcf_pipeline.py`)

```python
params = [
        (
            "test.vcf.gz",
            "fa758538aaea3a54de9dc302dd18e0d2"
        ),
        (
            "pindel_header.vcf",
            "9bd195a201d6b3317645ce5d44d40a2e"
        ),
        (
            "mutec2_indel.vcf.gz",
            "d03362c108c136ba45eb0d7de259a251"
        ),
        (
            "TCGA-DC-6682.vcf",
            "ea8daee9916c761edc9be74bc3eea475"
        )
]
```

### Updating Hash Values

When you intentionally change output format (bug fixes, new annotations), hashes will fail:

```bash
# Run tests with --update-vcf-hashes flag
pytest tests/integration/test_vcf_pipeline.py --update-vcf-hashes

# Output will show new hashes:
# New hash for test.vcf.gz: 3b0283fae9531d5874c95eab8cdcfe67
# New hash for mutec2_indel.vcf.gz: e97e468dcb509baf63710e341fdd8bdd
# ...

# Manually update the hashes in test_vcf_pipeline.py
# Then verify tests pass:
pytest tests/integration/test_vcf_pipeline.py
```

## Troubleshooting

### Common Issues

#### Issue: "TypeError: values expected to be 1-tuple, given len=2"

**Cause**: Trying to set INFO field with tuple/list when Number=1 in header.

**Fix**: Check INFO field definitions match your data:
```python
# If setting single value:
header.info.add('FIELD', 1, 'Integer', 'Description')
record.info['FIELD'] = 42  # Not (42,)

# If setting multiple values:
header.info.add('FIELD', '.', 'Integer', 'Description')
record.info['FIELD'] = (42, 43)  # Tuple/list OK
```

#### Issue: "OSError: unable to parse next record"

**Cause**: Malformed VCF file (duplicate #CHROM line, missing ##fileformat, etc.)

**Fix**: Always use pysam API to create VCF files in tests:
```python
# CORRECT
with pysam.VariantFile(temp_file.name, 'w', header=header) as vcf_out:
    vcf_out.write(record)

# WRONG - manual text manipulation
temp_file.write(str(header))
temp_file.write('#CHROM\tPOS\t...\n')
```

#### Issue: Hash validation fails after updating code

**Cause**: Output changed due to code modifications.

**Solution**:
1. Manually inspect `tests/output/` files to verify changes are correct
2. Run `pytest tests/integration/test_vcf_pipeline.py --update-vcf-hashes`
3. Update hash values in `test_vcf_pipeline.py`

### Known Failures (2 pre-existing issues)

#### 1. `test_sets_correct_position` (test_build_new_record.py:line 47)

**Issue**: BED coordinate system conversion inconsistency.

```python
# Test expects:
assert new_record.stop == 115  # BED END (1-based, non-inclusive)

# But gets:
assert new_record.stop == 114  # VCF stop (0-based)
```

**Status**: Pre-existing bug, not introduced by recent changes. 

#### 2. `test_repcn_for_heterozygous` (test_build_new_record.py:line 109)

**Issue**: Repeat count calculation for insertions.

```python
# Test case: A>AAT insertion in AT repeat
# Expected: REPCN should differ for heterozygous genotype (0/1)
# Actual: REPCN = (10, 10) - both alleles counted as 10 repeats

# Warning in logs:
# "Reference mismatch: VCF chr1:100 A>AAT, STR panel chr1:95-115 RU=AT"
```

**Status**: Pre-existing logic issue in how insertions are applied to repeat sequences.

## Continuous Integration

### GitHub Actions

The project uses GitHub Actions for automated testing on every push/PR:

```yaml
# .github/workflows/tests.yml
name: Tests
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.9'
      - run: pip install -r requirements_dev.txt
      - run: pytest tests/ --cov=src/strvcf_annotator
```

### Pre-commit Hooks

Install pre-commit hooks to run tests before committing:

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

## Test Coverage Summary

Current coverage by module (as of last test run):

| Module | Coverage | Lines | Missing |
|--------|----------|-------|---------|
| `api.py` | 100% | 67 | 0 |
| `cli.py` | 98% | 154 | 3 |
| `core/annotation.py` | 97% | 245 | 7 |
| `core/repeat_utils.py` | 99% | 187 | 2 |
| `core/str_reference.py` | 100% | 94 | 0 |
| `core/vcf_processor.py` | 96% | 203 | 8 |
| `parsers/base.py` | 100% | 23 | 0 |
| `parsers/generic.py` | 100% | 34 | 0 |
| `parsers/mutect.py` | 95% | 41 | 2 |
| `utils/validation.py` | 100% | 78 | 0 |
| **Overall** | **98%** | **1,126** | **22** |

**Note**: Coverage percentages are approximate. Run `pytest --cov` for exact current values.

## Contributing Tests

When adding new features or fixing bugs, follow this checklist:

- [ ] Write unit tests for new functions/methods
- [ ] Add integration test if feature affects end-to-end workflow
- [ ] Test edge cases (empty input, large files, invalid data)
- [ ] Test error conditions (missing files, malformed VCFs)
- [ ] Update hash values if output format changes
- [ ] Ensure all tests pass: `pytest tests/`
- [ ] Check coverage: `pytest tests/ --cov=src/strvcf_annotator`
- [ ] Add docstrings to test functions explaining what is tested
- [ ] Update this TESTS.md if adding new test categories

### Test Naming Conventions

- **Test files**: `test_<module_name>.py`
- **Test classes**: `class TestFunctionName:`
- **Test methods**: `def test_<specific_behavior>():`

Examples:
- ✅ `test_filters_variants_with_identical_genotypes()`
- ✅ `test_raises_value_error_on_missing_file()`
- ✅ `test_creates_correct_header_for_mutect2_vcf()`
- ❌ `test_function1()` (too vague)
- ❌ `test_stuff()` (meaningless)

---

**Last Updated**: Based on test suite as of 2025-11-04
**Test Statistics**: 156 passing / 158 total (98.7%)
**Maintainer**: See CONTRIBUTING.md for contact information
