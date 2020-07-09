from typing import Tuple, List
from unittest.mock import MagicMock

import pytest
from hypothesis import given
import hypothesis.strategies as st

from cigarco.cigar_utils import QUERY_CONSUMING_OPERATIONS, TARGET_CONSUMING_OPERATIONS
from cigarco.mapping import CMapper, Alignment
from test.test_cigar import decomposed_cigars


def test_mapper_invalid_initialization():
    with pytest.raises(TypeError):
        CMapper()


@pytest.fixture
def ex_alignment(scope="module"):
    return Alignment("tr1", "chr1", 3, "8M7D6M2I2M11D7M")


@pytest.fixture(scope="module")
def ex_cmapper():
    return CMapper(Alignment("tr1", "chr1", 3, "8M7D6M2I2M11D7M"))


@pytest.fixture(scope="module")
def ex_cmapper2():
    return CMapper(Alignment("tr2", "chr1", 10, "20M"))


def test_prefix_sum_array_lazy_initialization(ex_alignment):
    """
    This test checks that upon creation of the CMapper object the data-storing protected attributes _X_prefix_sums that store prefix sums for target and query consuming operation
        are not initialized

    Args:
        ex_alignment: a valid Alignment object
    """
    mapper = CMapper(alignment=ex_alignment)
    assert mapper._query_prefix_sums is None
    assert mapper._target_prefix_sums is None


def test_prefix_sum_array_lazy_computation(ex_alignment):
    """
    This test checks that upon creation of the CMapper object the data-storing protected attributes _X_prefix_sums that store prefix sums for target and query consuming operation
        are not initialized, but on first access to non-protected versions of said attributes, methods computing the prefix sum arrays are called, and the data is actually inferred

    Args:
        ex_alignment: a valid Alignment object
    """
    mapper = CMapper(alignment=ex_alignment)
    mapper.compute_target_prefix_sums = MagicMock(side_effect=mapper.compute_target_prefix_sums)
    mapper.compute_query_prefix_sums = MagicMock(side_effect=mapper.compute_query_prefix_sums)
    query_prefix_sums = mapper.query_prefix_sums
    target_prefix_sums = mapper.target_prefix_sums
    assert len(query_prefix_sums) > 0  # for the smallest CIGAR string (i.e., op count, operation) the length of the prefix sum array is going to be 2
    assert len(query_prefix_sums) == len(target_prefix_sums)  # must be the same for query and for target, as nonconsuming operations for each type are counted as 0
    mapper.compute_target_prefix_sums.assert_called()
    mapper.compute_query_prefix_sums.assert_called()
    assert query_prefix_sums is mapper._query_prefix_sums
    assert target_prefix_sums is mapper.target_prefix_sums
    # simple access to respective arrays shall not trigger subsequent computation method invocations
    mapper.query_prefix_sums
    mapper.query_prefix_sums
    mapper.target_prefix_sums
    assert mapper.compute_query_prefix_sums.call_count == 1
    assert mapper.compute_target_prefix_sums.call_count == 1


@given(data=st.lists(st.integers(min_value=0)))
def test_prefix_array_computation_logic(data):
    """
    inefficient O(n^2), though explicit computation of prefix sums in the test, but the algorithm implementation is more efficient O(n)
    """
    prefix_sums = CMapper.compute_prefix_sums(data)
    assert len(prefix_sums) == len(data)
    if len(prefix_sums) > 0:
        assert prefix_sums[-1] == sum(data)
    for i in range(len(data)):
        assert prefix_sums[i] == sum(data[:i + 1])


def test_prefix_sums_arrays_inference(ex_alignment):
    mapper = CMapper(ex_alignment)
    # the following prefix sums are based on the CIGAR string "8M7D6M2I2M11D7M" in the example_alignment object
    query_prefix_sums = [8, 8, 14, 16, 18, 18, 25]
    target_prefix_sums = [8, 15, 21, 21, 23, 34, 41]
    for true, inferred in zip(query_prefix_sums, mapper.query_prefix_sums):
        assert true == inferred
    for true, inferred in zip(target_prefix_sums, mapper.target_prefix_sums):
        assert true == inferred


@given(coord=st.integers(max_value=-1))
def test_coordinate_mapping_qt_invalid_negative(coord, ex_cmapper):
    """
    Only non-negative coordinate are allowed for transformation between query and target coordinate systems
    """
    with pytest.raises(ValueError):
        ex_cmapper.transform_coordinate(coord)


@given(coord=st.integers(min_value=25))
def test_coordinate_mapping_qt_invalid_over_query_length(coord, ex_cmapper):
    """
    Coordinates have to be <= the 0based length of the query sequence
    """
    with pytest.raises(ValueError):
        ex_cmapper.transform_coordinate(coord)


@given(coord=st.integers(min_value=0, max_value=24))
def test_coordinate_mapping_qt_ex1(coord, ex_cmapper):
    """
    Testing exactly the example 1 we have we have, where all the results are hand-checked and are comapred against
    """
    target_coord = ex_cmapper.transform_coordinate(coord)
    assert 3 <= target_coord <= 43
    if coord <= 7:
        assert target_coord == 3 + coord
    elif 7 < coord < 14:
        assert target_coord == 3 + coord + 7
    elif 14 <= coord < 16:
        assert target_coord == 23
    elif 16 <= coord < 18:
        assert target_coord in [24, 25]
    elif 18 <= coord < 25:
        assert target_coord == 3 + coord + 7 + 11 - 2


@given(coord=st.integers(min_value=0, max_value=19))
def test_coordinate_mapping_qt_ex2(coord, ex_cmapper2):
    target_coord = ex_cmapper2.transform_coordinate(coord)
    assert target_coord == 10 + coord


@st.composite
def alignment_tasks(draw):
    decomposed_cigar = draw(decomposed_cigars())
    data = list(zip(*decomposed_cigar))
    query_length = CMapper.compute_prefix_sums([cnt if op in QUERY_CONSUMING_OPERATIONS else 0 for cnt, op in data])[-1]
    coordinate = draw(st.integers(min_value=0, max_value=max(0, query_length - 1)))
    return decomposed_cigar, coordinate


@given(query_name=st.text(),
       target_name=st.text(),
       start=st.integers(min_value=0),
       alignment_task=alignment_tasks())
def test_coordinate_mapping_qt_random(query_name: str, target_name: str, start: int, alignment_task: Tuple[List[Tuple[int, str]], int]):
    """
    This test describes the overall logic between coordinate transformation from query to target coordinate systems.
    Note that we don't constrain the coordinate variable in the hypothesis generation strategy, as we are going to explicitly check for invalid cases
    """
    decomposed_cigar, coordinate = alignment_task
    decomposed_cigar = list(zip(*decomposed_cigar))
    cigar_string = "".join(f"{cnt}{op}" for cnt, op in decomposed_cigar)
    mapper = CMapper(Alignment(query_name, target_name, start, cigar_string))
    query_prefix_sums = CMapper.compute_prefix_sums([cnt if op in QUERY_CONSUMING_OPERATIONS else 0 for cnt, op in decomposed_cigar])
    target_prefix_sums = CMapper.compute_prefix_sums([cnt if op in TARGET_CONSUMING_OPERATIONS else 0 for cnt, op in decomposed_cigar])

    target_coordinate = mapper.transform_coordinate(coordinate)
    # general check that the query coordinate shall be transformed to the value of the target
    #   that does dont exceed the number of bases consumed in the target by a given alignment
    assert start <= target_coordinate <= start + target_prefix_sums[-1]
    true_target_coord = start
    last_target_only_operations = 0
    for cnt, op in decomposed_cigar:
        if coordinate < cnt and op in QUERY_CONSUMING_OPERATIONS:
            if op not in TARGET_CONSUMING_OPERATIONS:
                true_target_coord = max(true_target_coord - last_target_only_operations - 1, start)
            else:
                true_target_coord += coordinate
            assert true_target_coord == target_coordinate
            break
        else:
            if op in QUERY_CONSUMING_OPERATIONS:
                coordinate -= cnt
            if op in TARGET_CONSUMING_OPERATIONS:
                true_target_coord += cnt
            if cnt > 0 and op in TARGET_CONSUMING_OPERATIONS:
                if op not in QUERY_CONSUMING_OPERATIONS:
                    last_target_only_operations += cnt
                else:
                    last_target_only_operations = 0
    else:
        true_target_coord -= last_target_only_operations
        assert true_target_coord == target_coordinate


def test_mapping_qt_manual_insertion1():
    """
    Manual thought about test case where alignment ends with insertion, and simple deduction would lead to < start result
    """
    mapper = CMapper(Alignment("1", "1", 5, "2M10I"))
    assert mapper.transform_coordinate(8) == 6
    assert mapper.transform_coordinate(10) == 6


def test_mapping_qt_manual_insertion2():
    """
    Manual thought about test case where alignment ends with insertion, and simple deduction would lead to < start result
    """
    mapper = CMapper(Alignment("1", "1", 3, "2M2M10I"))
    assert mapper.transform_coordinate(8) == 6
    assert mapper.transform_coordinate(10) == 6


def test_mapping_qt_manual_insertion3():
    """
    Manual thought about test case where alignment ends with insertion, preceded by a deletion, and the transformation shall align with the left-most match w.r.t. insertion
    """
    mapper = CMapper(Alignment("1", "1", 3, "2M2M2D10I"))
    assert mapper.transform_coordinate(8) == 6
    assert mapper.transform_coordinate(10) == 6


def test_mapping_transform_coordinate_cache(ex_cmapper):
    ex_cmapper.transform_coordinate.cache_clear()
    ex_cmapper.transform_coordinate(5)
    assert ex_cmapper.transform_coordinate.cache_info().currsize == 1
    ex_cmapper.alignment = Alignment("1", "1", 3, "30M")
    assert ex_cmapper.transform_coordinate.cache_info().currsize == 0


def test_mapping_internal_data_structures_reset_on_alignment_attribute_update():
    mapper = CMapper(Alignment("1", "1", 3, "20M"))
    assert mapper.transform_coordinate(8) == 11
    assert mapper._query_prefix_sums is not None
    mapper.alignment = Alignment("1", "1", 3, "30M")
    assert mapper._query_prefix_sums is None


def test_mapping_transform_coordinate_cache_reset_on_alignment_attribute_update():
    mapper = CMapper(Alignment("1", "1", 3, "20M"))
    assert mapper.transform_coordinate(8) == 11
    assert mapper.transform_coordinate.cache_info().currsize == 1
    mapper.alignment = Alignment("1", "1", 3, "30M")
    assert mapper.transform_coordinate.cache_info().currsize == 0
