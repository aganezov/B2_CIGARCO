import pytest

from cigarco.mapping import Alignment


def test_empty_alignment_creation():
    with pytest.raises(TypeError):
        Alignment()


def test_basic_alignment_attributes():
    alignment = Alignment("tr1", "chr1", 0, "11M")
    assert alignment.query_name == "tr1"
    assert alignment.target_name == "chr1"
    assert alignment.start == 0
    assert alignment.cigar == "11M"


def test_invalid_start_in_alignment():
    with pytest.raises(ValueError):
        Alignment("tr1", "chr1", -1, "11M")

