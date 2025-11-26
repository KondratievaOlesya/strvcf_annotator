"""Utilities for repeat sequence operations."""

from typing import Dict, Union
import pandas as pd


def extract_repeat_sequence(str_row: Union[Dict, pd.Series]) -> str:
    """Reconstruct repeat sequence from STR metadata.

    Generates the full repeat sequence by repeating the repeat unit (RU)
    the calculated number of times (COUNT).

    Parameters
    ----------
    str_row : Dict or pd.Series
        STR region data containing 'RU' (repeat unit) and 'COUNT' (number of repeats)

    Returns
    -------
    str
        Full repeat sequence

    Examples
    --------
    >>> str_row = {'RU': 'CAG', 'COUNT': 5}
    >>> extract_repeat_sequence(str_row)
    'CAGCAGCAGCAGCAG'
    """
    return str_row["RU"] * int(str_row["COUNT"])


def count_repeat_units(sequence: str, motif: str) -> int:
    """Return the longest contiguous run of `motif` in `sequence`.

    The function looks for exact, non-overlapping copies of `motif` that occur
    consecutively and returns the maximum number of such copies in any run.

    This corresponds to how STR repeat counts are typically defined: the
    length of the longest perfect contiguous block of the repeat unit.

    Parameters
    ----------
    sequence : str
        DNA sequence to search.
    motif : str
        Repeat unit motif to count (e.g. 'A', 'CAG').

    Returns
    -------
    int
        Length of the longest contiguous run of `motif` in `sequence`.

    Raises
    ------
    ValueError
        If `motif` is empty or if either argument is not a string.

    Examples
    --------
    Perfect repeats
    ~~~~~~~~~~~~~~~
    >>> count_repeat_units("CAGCAGCAG", "CAG")
    3

    Imperfect tail
    ~~~~~~~~~~~~~~
    >>> count_repeat_units("CAGCAGCA", "CAG")
    2

    No repeats
    ~~~~~~~~~~
    >>> count_repeat_units("ATCG", "CAG")
    0

    Homopolymer runs
    ~~~~~~~~~~~~~~~~
    >>> count_repeat_units("ATAAAAA", "A")
    5
    >>> count_repeat_units("AAAATAAA", "A")
    4  # longest contiguous run is 'AAAA'

    Overlapping motifs
    ~~~~~~~~~~~~~~~~~~
    'AAAA' with motif 'AA' contains two non-overlapping copies: 'AA' 'AA'
    >>> count_repeat_units("AAAA", "AA")
    2
    """
    # Basic type and value validation
    if not isinstance(sequence, str):
        raise ValueError(f"'sequence' must be a string, got {type(sequence).__name__}")
    if not isinstance(motif, str):
        raise ValueError(f"'motif' must be a string, got {type(motif).__name__}")
    if motif == "":
        raise ValueError("'motif' must be a non-empty string")

    mlen = len(motif)
    slen = len(sequence)

    if mlen > slen or slen == 0:
        return 0

    max_run = 0
    i = 0

    while i <= slen - mlen:
        # Check if a run starts at position i
        if sequence[i : i + mlen] == motif:
            run_length = 0
            j = i
            # Count how many contiguous motif copies we have from position i
            while j <= slen - mlen and sequence[j : j + mlen] == motif:
                run_length += 1
                j += mlen

            if run_length > max_run:
                max_run = run_length

            # Skip past this entire run
            i = j
        else:
            i += 1

    return max_run

def apply_variant_to_repeat(
    pos: int, ref: str, alt: str, repeat_start: int, repeat_seq: str
) -> str:
    """Apply variant to repeat sequence, handling edge cases.

    Applies a variant (ref → alt) to the STR repeat sequence. Handles cases
    where the variant starts before the repeat region by clipping the variant
    appropriately.

    Parameters
    ----------
    pos : int
        Variant position (1-based VCF coordinate)
    ref : str
        Reference allele sequence
    alt : str
        Alternate allele sequence
    repeat_start : int
        Start position of repeat region (1-based)
    repeat_seq : str
        Full repeat sequence

    Returns
    -------
    str
        Mutated repeat sequence after applying the variant

    Examples
    --------
    >>> apply_variant_to_repeat(100, 'A', 'T', 100, 'AAAA')
    'TAAA'
    >>> apply_variant_to_repeat(100, 'AA', 'A', 100, 'AAAA')
    'AAA'

    Notes
    -----
    - Handles variants that start before the repeat region by clipping
    - Handles variants that extend beyond the repeat region
    - Preserves sequence integrity for complex indels
    """
    relative_pos = pos - repeat_start

    # Variant starts before the STR
    if relative_pos < 0:
        # Clip how much of the REF applies before the STR
        ref_clip = -relative_pos
        ref = ref[ref_clip:]  # Remove the part before STR
        alt = alt[ref_clip:]  # Same for ALT
        relative_pos = 0  # Mutation starts at beginning of repeat_seq

    # Apply mutation safely inside repeat_seq
    before_mut = repeat_seq[:relative_pos]
    after_mut = (
        repeat_seq[relative_pos + len(ref) :] if relative_pos + len(ref) <= len(repeat_seq) else ""
    )

    mutated = before_mut + alt + after_mut
    return mutated


def is_perfect_repeat(sequence: str, motif: str) -> bool:
    """Check if sequence is a perfect repeat of the motif.

    A perfect repeat means the sequence consists entirely of exact copies
    of the motif with no interruptions or variations.

    Parameters
    ----------
    sequence : str
        DNA sequence to check
    motif : str
        Repeat unit motif

    Returns
    -------
    bool
        True if sequence is a perfect repeat, False otherwise

    Examples
    --------
    >>> is_perfect_repeat('CAGCAGCAG', 'CAG')
    True
    >>> is_perfect_repeat('CAGCAGCA', 'CAG')
    False
    """
    if not sequence or not motif:
        return False

    count = count_repeat_units(sequence, motif)
    return sequence == motif * count
