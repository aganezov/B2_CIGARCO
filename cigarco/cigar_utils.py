ALLOWED_OPERATIONS = {'M', "I", "D", "N", "S", "H", "P", "=", "X"}
QUERY_CONSUMING_OPERATIONS = ALLOWED_OPERATIONS - {"D", "N", "H", "P"}
TARGET_CONSUMING_OPERATIONS = ALLOWED_OPERATIONS - {"I", "S", "H", "P"}


def is_valid_cigar(cigar: str) -> bool:
    """
    A CIGAR string is valid if it is non-empty string with alternating numbers and supporting single character operations

    Args:
        cigar (str): a CIGAR encoded (doc: https://samtools.github.io/hts-specs/SAMv1.pdf , page 7) alignment string

    Returns:
        valid (bool): a flag indicating if a cigar string is valid one.
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
