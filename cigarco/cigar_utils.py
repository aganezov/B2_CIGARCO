import re
from typing import Tuple, List

ALLOWED_OPERATIONS = {'M', "I", "D", "N", "S", "H", "P", "=", "X"}
QUERY_CONSUMING_OPERATIONS = ALLOWED_OPERATIONS - {"D", "N", "H", "P"}
TARGET_CONSUMING_OPERATIONS = ALLOWED_OPERATIONS - {"I", "S", "H", "P"}

# A regex designed to capture in a named way pairs of number encoded counts for operations, as well as the operations themselves
#   meant to be used in a `re.finditer` fashion
# This regex is going to be reused constantly, so its better to compute it once on the module load
CIGAR_REGEX = re.compile(r'((?P<count>\d+)(?P<operation>[MIDNSHP=X]))')


def is_valid_cigar(cigar: str) -> bool:
    """
    A CIGAR string is valid if it is non-empty string with alternating numbers and supporting single character operations
    Works in a single pass with O(n) complexity, where n is the size of the input string; requires O(1) extra memory.

    Args:
        cigar (str): a CIGAR encoded (doc: https://samtools.github.io/hts-specs/SAMv1.pdf , page 7) alignment string

    Returns:
        valid (bool): a flag indicating if a cigar string is valid one.

    Examples:
        >>> is_valid_cigar("M")
        False

        >>> is_valid_cigar("10")
        False

        >>> is_valid_cigar("10M10")
        False

        >>> is_valid_sigar("10M10D")
        True
    """
    # can't be an empty one
    if not cigar:
        return False
    prev = None
    for char in cigar:

        # Every entry has to be either a digit or a supporter single-char encoded operation
        if not (char.isdigit() or char in ALLOWED_OPERATIONS):
            return False

        # CIGAR string must start with a digit
        if char in ALLOWED_OPERATIONS and (not prev or not prev.isdigit()):
            return False

        prev = char

    # CIGAR string must end with an operation
    if prev not in ALLOWED_OPERATIONS:
        return False

    return True


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """
    Parses a given CIGAR encoded string into a list of tuples (count, operation)
    No internal checks for the validity of the input CIGAR string are made

    Args:
        cigar (str): CIGAR encoded string

    Returns:
        a list of tuples CIGAR operations (str) and their counts (List[Tuple[int, str]])

    Examples:
        >>> parse_cigar("10M11D")
        [(10, "M"), (11, "D")]

        >>> parse_cigar("1M115I10M")
        [(1, "M"), (115 "I"), (10, "M")]
    """
    result: List[Tuple[int, str]] = []
    for entry in re.finditer(CIGAR_REGEX, cigar):
        result.append((int(entry.groupdict()["count"]), entry.groupdict()["operation"]))
    return result
