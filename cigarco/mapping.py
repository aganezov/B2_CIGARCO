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


@dataclass
class CMapper(object):
    """
    The main object that provides coordinate conversion for a given Alignment object.
    The main computational idea behind the CMapper machinery is the prefix sum arrays, which are used to efficiently convert the query coordinate to the target coordinate

    Several optimization techniques are utilized to speed up the process:
        * prefix sum arrays are computed only once on first access via a combination of private attributes and property access control
        * alignment object reassignment (controlled via property setter) invalidates the computed prefix sums ensuring that coordinate mapping is always done
            w.r.t. current alignment; Note: changes in alignment object are nto tracked, though this can be implemented in the future via a descriptor-protocol on alignment
        * cashing is utilized for mapping function

    Overall complexity of the CMapper is dependent on the size n of the CIGAR alignment string (where n refers to the number of operations in a parsed CIGAR string):
        * prefix sum construction takes O(n); on first access only
        * coordinate mapping takes O(log(n)) steps at worst (and O(1) for repeated lookups)

    Args:
        alignment (Alignments): an alignment object for which the coordinate conversion (i.e., mapping) is performed. Alignment object is assumed to be valid
    """
    alignment: Alignment
