from dataclasses import dataclass


@dataclass
class Alignment(object):
    query_name: str
    target_name: str
    start: int
    cigar: str
