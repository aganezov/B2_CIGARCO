from copy import deepcopy
from typing import List

from hypothesis import given

from cigarco.mapping import CManager, Alignment, CMapper
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
