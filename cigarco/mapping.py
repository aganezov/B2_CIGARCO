from dataclasses import dataclass
from typing import List, Tuple

from cigarco.cigar_utils import is_valid_cigar, parse_cigar, TARGET_CONSUMING_OPERATIONS, QUERY_CONSUMING_OPERATIONS


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
    _query_prefix_sums: List[int] = None
    _target_prefix_sums: List[int] = None

    @property
    def query_prefix_sums(self) -> List[int]:
        """ Implements a descriptor-like protocol for accessing prefix sums arrays for query consuming operations in CIGAR string
            computation is only invoked if the underlying data-holding attribute is set to None, which can happen either before the first access to the prefix sums array,
            or if the alignment object has been changed and prefix sums arrays have been invalidated

        Returns:
            prefix sums (List[int]) for all operations in alignment CIGAR string
        """
        if self._query_prefix_sums is None:
            self.compute_query_prefix_sums()
        return self._query_prefix_sums

    def compute_query_prefix_sums(self):
        """
        Computes prefix sums array for query consuming operations in the alignment object's CIGAR string

        TODO: rewrite with iterators approach
        """
        cigar_operations: List[Tuple[int, str]] = parse_cigar(self.alignment.cigar)
        query_op_cnts = [cnt if op in QUERY_CONSUMING_OPERATIONS else 0 for cnt, op in cigar_operations]
        self._query_prefix_sums = self.compute_prefix_sums(query_op_cnts)

    @property
    def target_prefix_sums(self) -> List[int]:
        """ Implements a descriptor-like protocol for accessing prefix sums arrays for target consuming operations in CIGAR string
            computation is only invoked if the underlying data-holding attribute is set to None, which can happen either before the first access to the prefix sums array,
            or if the alignment object has been changed and prefix sums arrays have been invalidated

        Returns:
            prefix sums (List[int]) for all operations in alignment CIGAR string
        """
        if self._target_prefix_sums is None:
            self.compute_target_prefix_sums()
        return self._target_prefix_sums

    def compute_target_prefix_sums(self):
        """
        Computes prefix sums array for target consuming operations in the alignment object's CIGAR string

        TODO: rewrite with iterators approach
        """
        cigar_operations: List[Tuple[int, str]] = parse_cigar(self.alignment.cigar)
        target_op_cnts = [cnt if op in TARGET_CONSUMING_OPERATIONS else 0 for cnt, op in cigar_operations]
        self._target_prefix_sums = self.compute_prefix_sums(target_op_cnts)   # TODO: remove dummy implementation

    @staticmethod
    def compute_prefix_sums(values: List[int]) -> List[int]:
        """
        Stateless (thus staticmethod) utility function that computes prefix sums array for a given integer array.
        For a given array A a prefix sum array A' is defined as follows:
            A'[0] = 0
            for i > 0: A'[i] = sum(A[0] ... A[i-1])
        Provided implementation works in linear O(n) time, where `n` is the length of the input array

        Args:
            values (List[int]): a list of integers

        Returns:
            prefix sums (List[int]): a list a prefix sums for the input list

        Examples:
            >>> CMapper.compute_prefix_sums([1,2,3])
            [0,1,3,6]

            >>> CMapper.compute_prefix_sums([])
            [0]
        """
        result = [0]
        for v in values:
            result.append(result[-1] + v)
        return result

    def transform_coordinate(self, source_coordinate: int) -> int:
        """
        The main method to be invoked for coordinate transformation from query -> target coordinate systems.
        Uses precomputed (or lazily evaluated on the first invocation) prefix sums arrays for query/target consuming operations and binary search for efficient lookup.
        Also utilizes memoization to reduce computational load in case of identical transformation requests (cache is invalidated if the alignment object is altered)

        Coordinates that map to insertion in query (w.r.t. target) are transformed to the coordinate of the insertion start
            (i.e., the last target coordinate before the insertion seq)

        Args:
            source_coordinate (int): source coordinate in query coordinate system to be translated to the target coordinate system

        Returns:
            target coordinate (int): a position that the argument source coordinate maps to in the target sequence

        Raises:
            ValueError: if the input source coordinate is negative or greater than the length of the query sequence
        """
        if source_coordinate < 0 or source_coordinate > self.query_prefix_sums[-1] - 1:
            # last value in prefix sums array is the length of the query, but we need to account for the 0-based index we
            raise ValueError(f"Can't transform coordinate {source_coordinate}, outside of query coordinate system")
