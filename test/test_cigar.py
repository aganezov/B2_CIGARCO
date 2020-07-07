import pytest

from cigarco.mapping import Alignment


def test_empty_alignment_creation():
    with pytest.raises(TypeError):
        Alignment()
