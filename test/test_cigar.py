from hypothesis import given, settings
import hypothesis.strategies as st
import string

from cigarco.cigar_utils import is_valid_cigar, ALLOWED_OPERATIONS, parse_cigar


def test_empty_cigar_validity():
    assert not is_valid_cigar("")


@given(number=st.integers(min_value=0), operation=st.text(alphabet=ALLOWED_OPERATIONS, min_size=2, max_size=100))
def test_invalid_cigar_consecutive_letters(number, operation):
    assert not is_valid_cigar(f"{number}{operation}")


@given(number=st.integers(min_value=0))
def test_invalid_cigar_no_letters(number):
    assert not is_valid_cigar(f"{number}")


@given(operation=st.text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1))
def test_invalid_cigar_no_starting_number(operation):
    assert not is_valid_cigar(f"{operation}")


@given(n1=st.integers(min_value=0),
       operation=st.text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1),
       n2=st.integers(min_value=0))
def test_invalid_cigar_no_ending_operation(n1, operation, n2):
    assert not is_valid_cigar(f"{n1}{operation}{n2}")


@given(n1=st.integers(min_value=0),
       operation=st.text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1),
       n2=st.integers(min_value=0),
       unsupported_operation=st.text(alphabet=set(string.ascii_uppercase) - ALLOWED_OPERATIONS, min_size=1, max_size=1))
def test_invalid_cigar_unsupported_operation_character(n1, operation, n2, unsupported_operation):
    assert not is_valid_cigar(f"{n1}{operation}{n2}{unsupported_operation}")


@st.composite
def cigar_decomposition(draw):
    counts = draw(st.lists(elements=st.integers(min_value=0), min_size=1, max_size=1000000000))
    operations = draw(st.lists(elements=st.text(alphabet=ALLOWED_OPERATIONS, min_size=1, max_size=1), min_size=1, max_size=1000000000))
    size = min(len(counts), len(operations))
    return counts[:size], operations[:size]


@given(values=cigar_decomposition())
def test_cigar_split(values):
    counts, operations = values
    assert len(counts) == len(operations)
    cigar_string_data = [(c, o) for c, o in zip(counts, operations)]
    cigar_string = "".join(f"{c}{o}" for c, o in cigar_string_data)
    parsed_data = parse_cigar(cigar_string)
    assert len(cigar_string_data) == len(parsed_data)
    for true, inferred in zip(cigar_string_data, parsed_data):
        assert true[0] == inferred[0]
        assert true[1] == inferred[1]
