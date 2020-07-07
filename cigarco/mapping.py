from dataclasses import dataclass

from cigarco.cigar_utils import is_valid_cigar


@dataclass
class Alignment(object):
    """
    An alignment class an instance of which contains all the information about a given query->target alignment

    Args:
        query_name (str): a name of the query which alignment we represent
        target_name (str): a name of the target for the represented alignment
        start (int): a start coordinate for the query alignment w.r.t. to the target
        cigar (str): a CIGAR encoded (doc: https://samtools.github.io/hts-specs/SAMv1.pdf , page 7) alignment of the query to the target

    Example:
        Alignment("tr1", "chr1" 0, "11M")
    """
    query_name: str
    target_name: str
    start: int          # no default value of 0 specified as different schools of thought may have 0 or 1 as defaults, and explicit is better than implicit
    cigar: str

    def __post_init__(self):
        """
        Ensuring the validity of the created alignment abject, as subsequent methods that accept the alignment object expect a valid alignment
        Raises:
            ValueError: if the start coordinate is negative
        """
        if self.start < 0:
            raise ValueError(f"incorrect start coordinate {self.start}. Must be a non-negative integer")
        if not is_valid_cigar(self.cigar):
            raise ValueError(f"invalid CIGAR string '{self.cigar}'")
