import pytest
from hypothesis import given
import hypothesis.strategies as st

from cigarco.mapping import Alignment


def test_empty_alignment_creation():
    with pytest.raises(TypeError):
        Alignment()


@given(query=st.text(), target=st.text(), start=st.integers(min_value=0))
def test_basic_alignment_attributes(query, target, start):
    cigar = "10M"   # manual valid cigar to ensure we pass the validity check
    alignment = Alignment(query, target, start, cigar)
    assert alignment.query_name == query
    assert alignment.target_name == target
    assert alignment.start == start
    assert alignment.cigar == cigar


@given(start=st.integers(max_value=-1))
def test_invalid_start_in_alignment(start):
    with pytest.raises(ValueError):
        Alignment("tr1", "chr1", start, "11M")


def test_invalid_cigar_argument():
    with pytest.raises(ValueError):
        Alignment("tr1", "chr1", 0, "10")

