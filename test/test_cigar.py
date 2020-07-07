import pytest

from cigarco.cigar_utils import is_valid_cigar


@pytest.mark.xfail(reason="empty CIGAR string is not valid")
def test_empty_cigar_validity():
    assert is_valid_cigar("")


@pytest.mark.xfail(reason="CIGAR string does nto allow for consecutive letters")
def test_invalid_cigar_consecutive_letters():
    assert is_valid_cigar("11MM")


@pytest.mark.xfail(reason="CIGAR string does not allow for consecutive letters")
def test_invalid_cigar_no_letters():
    assert is_valid_cigar("11")


@pytest.mark.xfail(reason="CIGAR string must have a number encoding the amount of time a subsequently specified operation is encoded")
def test_invalid_cigar_no_starting_number():
    assert is_valid_cigar("M")


@pytest.mark.xfail(reason="CIGAR string must have an operation encoded  as a single letter after every provided number")
def test_invalid_cigar_no_ending_operation():
    assert is_valid_cigar("11M10")


@pytest.mark.xfail(reason="CIGAR string only supports MIDNSHP=X operation encoding characters")
def test_invalid_cigar_unsupported_operation_character():
    assert is_valid_cigar("10A11M")



