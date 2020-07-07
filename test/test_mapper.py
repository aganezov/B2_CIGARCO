import pytest

from cigarco.mapping import CMapper


def test_mapper_invalid_initialization():
    with pytest.raises(TypeError):
        CMapper()


