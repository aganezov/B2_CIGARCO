from unittest.mock import MagicMock

import pytest

from cigarco.mapping import CMapper, Alignment


def test_mapper_invalid_initialization():
    with pytest.raises(TypeError):
        CMapper()


@pytest.fixture
def ex_alignment():
    return Alignment("tr1", "chr1", 3, "8M7D6M2I2M11D7M")


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
    assert len(query_prefix_sums) > 1           # for the smallest CIGAR string (i.e., op count, operation) the length of the prefix sum array is going to be 2
    assert len(query_prefix_sums) == len(target_prefix_sums)    # must be the same for query and for target, as nonconsuming operations for each type are counted as 0
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


