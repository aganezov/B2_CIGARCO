import bisect
from collections import namedtuple
from dataclasses import dataclass, field
from functools import lru_cache
from typing import List, Tuple, Set, Optional, Dict
from itertools import accumulate

import typing

from cigarco.cigar_utils import is_valid_cigar, parse_cigar, TARGET_CONSUMING_OPERATIONS, QUERY_CONSUMING_OPERATIONS


@dataclass(frozen=True, eq=True)
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
    start: int  # no default value of 0 specified as different schools of thought may have 0 or 1 as defaults, and explicit is better than implicit
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


class AlignmentDescriptor(object):
    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner=None) -> Optional[Alignment]:
        return instance.__dict__.get(self.name)

    def __set__(self, instance, value):
        if isinstance(instance, CMapper):
            instance._query_prefix_sums = None
            instance._target_prefix_sums = None
            instance._matching_backtracking = None
            instance.transform_coordinate.cache_clear()
        instance.__dict__[self.name] = value


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
    alignment: Alignment = AlignmentDescriptor()
    _alignment: Alignment = field(init=False, repr=False, compare=False)
    _query_prefix_sums: Optional[List[int]] = field(init=False, repr=False, compare=False)
    _target_prefix_sums: Optional[List[int]] = field(init=False, repr=False, compare=False)
    _matching_backtracking: Optional[List[int]] = field(init=False, repr=False, compare=False)

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

    @property
    def matching_backtracking(self) -> List[int]:
        if self._matching_backtracking is None:
            self._matching_backtracking = self.compute_matching_backtracking()
        return self._matching_backtracking

    def compute_matching_backtracking(self) -> List[int]:
        last_match_index: int = -1
        result: List[int] = []
        qt_consuming_operations: Set[str] = QUERY_CONSUMING_OPERATIONS & TARGET_CONSUMING_OPERATIONS
        for index, (cnt, op) in enumerate(parse_cigar(self.alignment.cigar)):
            if cnt > 0 and op in qt_consuming_operations:
                last_match_index = index
            result.append(last_match_index)
        return result

    def compute_target_prefix_sums(self):
        """
        Computes prefix sums array for target consuming operations in the alignment object's CIGAR string

        TODO: rewrite with iterators approach
        """
        cigar_operations: List[Tuple[int, str]] = parse_cigar(self.alignment.cigar)
        target_op_cnts = [cnt if op in TARGET_CONSUMING_OPERATIONS else 0 for cnt, op in cigar_operations]
        self._target_prefix_sums = self.compute_prefix_sums(target_op_cnts)  # TODO: remove dummy implementation

    @staticmethod
    def compute_prefix_sums(values: List[int]) -> List[int]:
        """
        Stateless (thus staticmethod) utility function that computes prefix sums array for a given integer array.
        For a given array A a prefix sum array A' is defined as follows:
            for i >= 0: A'[i] = sum(A[0] ... A[i])
        Provided implementation works in linear O(n) time, where `n` is the length of the input array

        Args:
            values (List[int]): a list of integers

        Returns:
            prefix sums (List[int]): a list a prefix sums for the input list

        Examples:
            >>> CMapper.compute_prefix_sums([1,2,3])
            [1,3,6]

            >>> CMapper.compute_prefix_sums([])
            []
        """
        return list(accumulate(values))

    @lru_cache(maxsize=None)
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

        Examples:
            >>> CMapper(Alignment("tr1", "chr1", 3, "8M7D6M2I2M11D7M")).transform_coordinate(4)
            7

            >>> CMapper(Alignment("tr1", "chr1", 3, "8M7D6M2I2M11D7M")).transform_coordinate(13)
            23

            >>> CMapper(Alignment("tr2", "chr2", 10, "20M")).transform_coordinate(0)
            10

            >>> CMapper(Alignment("tr2", "chr2", 10, "20M")).transform_coordinate(10)
            20
        """
        if source_coordinate < 0 or source_coordinate > max(0, self.query_prefix_sums[-1] - 1):
            # last value in prefix sums array is the length of the query, but we need to account for the 0-based index
            raise ValueError(f"Can't transform coordinate {source_coordinate}, outside of query coordinate system")
        # identifies the last operation that would have consumed x <= source_coordinate bases in th query
        operation_index: int = bisect.bisect_right(self.query_prefix_sums, source_coordinate)

        # special edge case where the alignment cigar string has no query consuming operations, in which case we default to the beginning of the alignment
        # while highly improbable -- still allowed by the CIGAR specification in the SAM format
        if source_coordinate == 0 and self.query_prefix_sums[-1] == 0:
            return self.alignment.start

        # sanity check, does not waste resources at all, but if really needed, can be avoided with with -O flag in execution
        assert 0 <= operation_index <= len(self.query_prefix_sums) - 1

        # we get the operation index for the last stretch where there was a match between query an alignment,
        #   this is needed for ensuring that coordinates in non-target-consuming operations (i.e., insertions) map to the left-closest matching position
        last_matching_index: int = self.matching_backtracking[operation_index]

        # computing how much query has been consumed by the latest operation not covering the queried coordinate,
        #   this is required for figure out how much of a non-consumed query we are left with
        query_consumed: int = 0 if operation_index == 0 else self.query_prefix_sums[operation_index - 1]

        # computing the amount of target-matching query we are left with
        #   of the last_matching index and operation indexes don't match, we are in the non-target-consuming are need to set remaining query length to -1,
        #   to ensure left-padded insertion coordinate mapping
        query_remaining = -1 if operation_index != last_matching_index else (source_coordinate - query_consumed)

        # if we are in a matching operation, we need to decrement the matching operation index by 1 to ensure that we correctly calculate consumed target sequence
        #   (up until the identified operation)
        last_matching_index -= int(last_matching_index == operation_index)

        # target is only consumed to the last matching operation (not counting the one in which the query coordinate lies)
        target_consumed: int = 0 if last_matching_index < 0 else self.target_prefix_sums[last_matching_index]

        # we need to ensure that we don't end up with negative offset, which can come from weird valid CIGAR strings and the setup above (e.g., "I5")
        result: int = self.alignment.start + max(target_consumed + query_remaining, 0)
        return result

    def __hash__(self):
        return hash(self.alignment)


@dataclass(frozen=True, eq=True)
class TransformedCoordinate(object):
    """
    Simple holder class for the result of coordinate transformation where both the target seq and the transformed coordinate values are held
    Allows for nice data ncapsulation, while not eating almost any extra space, because of the __slots__ usage
    """
    __slots__ = ('seq_name', 'coordinate')
    seq_name: str
    coordinate: int


@dataclass
class CManager(object):
    """ Main class that manages storage of multiple alignments (one alignment per query) and allows for efficient coordinate transformation from query coordinate system
            to that of the alignment target one
        Outsources the actual computations to the CMapper class, that is created for every added alignment
        Keeps track of added alignments so as to when an identical alignment is added no recomputation is going to be performed

        First query for a given alignment takes O(m) + O(log(m)), where m is the number of operation in the alignment CIGAR string
        Subsequent queries of previously UNqueried values take O(log(m))
        Queries are cached, so subsequent queries of previously queried values take O(1)

        Examples:
            # setup
             >>> m = CManager()
             >>> m.add_alignment(Alignment("TR1", "CHR1", 3, "8M7D6M2I2M11D7M"))
             >>> m.add_alignment(Alignment("TR2", "CHR2", 10, "20M"))
             # quering
             >>> m.transform_coordinate("TR1", 4)
             TransformedCoordinate(seq_name="CHR1", coordinate=7)
             >>> m.transform_coordinate("TR2", 0)
             TransformedCoordinate(seq_name="CHR2", coordinate=10)

    """
    alignments_by_query_ids: Dict[str, CMapper] = field(default_factory=lambda: {})

    def add_alignment(self, alignment: Alignment):
        """
        Wrapper method that adds a new Alignment instance to the internal structure of the Manager, and wraps the alignment object into CMapper object
        When addition of a duplicate alignment is attempted, no action is performed, thus keeping the potentially computed coordinate transformation data structures intact

        Args:
            alignment (Alignment): an instance of an Alignment class
        """
        if alignment.query_name in self.alignments_by_query_ids:
            mapper = self.alignments_by_query_ids[alignment.query_name]
            if mapper.alignment != alignment:
                self.alignments_by_query_ids[alignment.query_name].alignment = alignment
        else:
            self.alignments_by_query_ids[alignment.query_name] = CMapper(alignment)

    def transform_coordinate(self, query_name: str, source_coordinate: int) -> TransformedCoordinate:
        """ The main computational method for coordinate transformation for a given position in a specified query
        On its own just checks if the mapper for a given query exists and then outsources the actual new coordinate value computation to it
        Caching is implemented at the level of a mapper, and not here, as with addition of (new) alignments only parts of cache would ideally be invalidated (for respective query),
            but this is not feasible with basic lru_cache, though at the level of the mapper this is exactly what is implemented

        Args:
            query_name (str): name of the query for which alignment the coordinate transformation is going to take place
            source_coordinate (int): a coordinate on specified query, which is going to be transformed into a alignment target corodinate system

        Returns:
            transformed coordinate (TransformedCoordinate): a dataclass (akin to pair tuple) of a name of a target sequence,
                and a transformed query coordinate in target coordinate system

        Raises:
            ValueError: if the supplied query does not have an alignment, or if the coordinate transformation fails (propagated from CMapper.transform_coordinate method)
        """
        if query_name not in self.alignments_by_query_ids:
            raise ValueError(f"Attempted to transform coordinate {source_coordinate} for query '{query_name}', but no alignments for '{query_name}' exist")
        mapper: CMapper = self.alignments_by_query_ids[query_name]
        target_seq_name: str = mapper.alignment.target_name
        result_coordinate: int = mapper.transform_coordinate(source_coordinate)
        return TransformedCoordinate(target_seq_name, result_coordinate)
