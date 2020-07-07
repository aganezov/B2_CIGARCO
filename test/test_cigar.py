from hypothesis import given
from hypothesis.strategies import text, integers
import string

from cigarco.cigar_utils import is_valid_cigar, ALLOWED_OPERATIONS


def test_empty_cigar_validity():
    assert not is_valid_cigar("")


@given(number=integers(min_value=0), operation=text(alphabet=ALLOWED_OPERATIONS, min_size=2, max_size=100))
def test_invalid_cigar_consecutive_letters(number, operation):
    assert not is_valid_cigar(f"{number}{operation}")


@given(number=integers(min_value=0))
def test_invalid_cigar_no_letters(number):
    assert not is_valid_cigar(f"{number}")


@given(operation=text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1))
def test_invalid_cigar_no_starting_number(operation):
    assert not is_valid_cigar(f"{operation}")


@given(n1=integers(min_value=0),
       operation=text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1),
       n2=integers(min_value=0))
def test_invalid_cigar_no_ending_operation(n1, operation, n2):
    assert not is_valid_cigar(f"{n1}{operation}{n2}")


@given(n1=integers(min_value=0),
       operation=text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1),
       n2=integers(min_value=0),
       unsupported_operation=text(alphabet=set(string.ascii_uppercase) - ALLOWED_OPERATIONS, min_size=1, max_size=1))
def test_invalid_cigar_unsupported_operation_character(n1, operation, n2, unsupported_operation):
    assert not is_valid_cigar(f"{n1}{operation}{n2}{unsupported_operation}")
