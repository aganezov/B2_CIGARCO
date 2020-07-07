import pytest

from cigarco.cigar_utils import is_valid_cigar


@pytest.mark.xfail(reason="empty CIGAR string is not valid")
def test_empty_cigar_validity():
    assert is_valid_cigar("")


@pytest.mark.xfail(reason="CIGAR string does nto allow for consecutive letters")
def test_invalid_cigar_consecutive_letters():
    assert is_valid_cigar("11MM")


@pytest.mark.xfail(reason="CIGAR string does nto allow for consecutive letters")
def test_invalid_cigar_no_letters():
    assert is_valid_cigar("11")

