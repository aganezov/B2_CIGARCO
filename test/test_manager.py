from collections import namedtuple
from copy import deepcopy
from typing import List, Tuple

import pytest
import typing
from hypothesis import given

from cigarco.cigar_utils import QUERY_CONSUMING_OPERATIONS
from cigarco.mapping import CManager, Alignment, CMapper, TransformedCoordinate
import hypothesis.strategies as st

from test.test_cigar import decomposed_cigars


def test_manager_creation():
    manager = CManager()


@st.composite
def alignments(draw):
    decomposed_cigar = draw(decomposed_cigars())
    cigar_string = "".join(f"{cnt}{op}" for cnt, op in zip(*decomposed_cigar))
    query_name = draw(st.text())
    target_name = draw(st.text())
    start = draw(st.integers(0))
    return Alignment(query_name, target_name, start, cigar_string)


@given(ex_alignments=st.lists(alignments()))
def test_alignment_addition(ex_alignments: List[Alignment]):
    manager: CManager = CManager()
    unique_query_alignments = {a.query_name: a for a in ex_alignments}
    for alignment in unique_query_alignments.values():
        manager.add_alignment(alignment)
    assert len(set(manager.alignments_by_query_ids) & set(unique_query_alignments)) == len(manager.alignments_by_query_ids)
    for k, v in manager.alignments_by_query_ids.items():
        assert isinstance(k, str)
        assert isinstance(v, CMapper)


@given(alignment=alignments())
def test_addition_of_existing_alignment(alignment: Alignment):
    manager: CManager = CManager()
    manager.add_alignment(alignment)
    # not the right way to do this, only done for testing
    manager.alignments_by_query_ids[alignment.query_name].transform_coordinate(0)
    assert manager.alignments_by_query_ids[alignment.query_name]._query_prefix_sums is not None
    mapper = manager.alignments_by_query_ids[alignment.query_name]
    manager.add_alignment(deepcopy(alignment))
    assert manager.alignments_by_query_ids[alignment.query_name]._query_prefix_sums is not None
    assert mapper is manager.alignments_by_query_ids[alignment.query_name]


@given(alignment=alignments())
def test_preservation_of_mapper(alignment: Alignment):
    manager: CManager = CManager()
    manager.add_alignment(alignment)
    mapper = manager.alignments_by_query_ids[alignment.query_name]
    new_alignment = Alignment(alignment.query_name, alignment.target_name, alignment.start + 1, alignment.cigar)
    manager.add_alignment(new_alignment)
    assert manager.alignments_by_query_ids[alignment.query_name] is mapper


@st.composite
def transformation_tasks(draw):
    decomposed_cigar = draw(decomposed_cigars())
    cigar_string = "".join(f"{cnt}{op}" for cnt, op in zip(*decomposed_cigar))
    query_length = 0
    for cnt, op in zip(*decomposed_cigar):
        if op in QUERY_CONSUMING_OPERATIONS:
            query_length += cnt
    query_name = draw(st.text())
    target_name = draw(st.text())
    start = draw(st.integers(0))
    coordinate = draw(st.integers(min_value=0, max_value=max(0, query_length - 1)))
    return Alignment(query_name, target_name, start, cigar_string), coordinate


@given(transformation_task=transformation_tasks())
def test_coordinate_transformation_qt(transformation_task: Tuple[Alignment, int]):
    manager: CManager = CManager()
    alignment, coordinate = transformation_task
    manager.add_alignment(alignment)
    transformed_coordinate: TransformedCoordinate = manager.transform_coordinate(alignment.query_name, coordinate)
    assert transformed_coordinate.seq_name == alignment.target_name
    assert transformed_coordinate.coordinate == manager.alignments_by_query_ids[alignment.query_name].transform_coordinate(coordinate)


@given(transformation_task=transformation_tasks())
def test_coordinate_transformation_qt_invalid_query_name(transformation_task: Tuple[Alignment, int]):
    manager: CManager = CManager()
    alignment, coordinate = transformation_task
    manager.add_alignment(alignment)
    with pytest.raises(ValueError):
        manager.transform_coordinate(alignment.query_name + "t", coordinate)


@given(transformation_task=transformation_tasks())
def test_coordinate_transformation_qt_invalid_coordinate_negative(transformation_task: Tuple[Alignment, int]):
    manager: CManager = CManager()
    alignment, coordinate = transformation_task
    manager.add_alignment(alignment)
    with pytest.raises(ValueError):
        manager.transform_coordinate(alignment.query_name, -1)


@given(transformation_task=transformation_tasks())
def test_coordinate_transformation_qt_invalid_coordinate_greater_than_query_length(transformation_task: Tuple[Alignment, int]):
    manager: CManager = CManager()
    alignment, coordinate = transformation_task
    manager.add_alignment(alignment)
    # cheating a bit, but ok for purposes of the test
    query_length = manager.alignments_by_query_ids[alignment.query_name].query_prefix_sums[-1]
    with pytest.raises(ValueError):
        manager.transform_coordinate(alignment.query_name, max(query_length, 1))
