from dataclasses import dataclass


@dataclass
class Alignment(object):
    query_name: str
    target_name: str
    start: int
    cigar: str

    def __post_init__(self):
        if self.start < 0:
            raise ValueError(f"incorrect start coordinate {self.start}. Must be a non-negative integer")
