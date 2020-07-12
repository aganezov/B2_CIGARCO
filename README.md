# CIGAR Coordinate (re)mapper

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python package](https://github.com/aganezov/B2_CIGARCO/workflows/Python%20package/badge.svg)](https://github.com/aganezov/B2_CIGARCO/actions)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/)
 

A toolkit for CLI- and API-based (re)mapping of alignment coordinates to/from reference coordinates based on CIGAR-encoded alignments.

### Contents:
* [Design and complexity](#design-and-complexity)
* [Installation](#installation)
* [Usage](#usage)
* [Future directions](#future-directions)

### Design and implementation

The main computational task the `CIGARCO` solves is the computationally efficient transformation of coordinate positions from one coordinate system (e.g., query) into another one (e.g., target).
The relationship between query and target sequences is determined by an alignment in a form of a CIGAR string (more details about CIGAR: [sam format](https://samtools.github.io/hts-specs/SAMv1.pdf), page 7).
`CIGARCO` can be both used as a library for integration into existing projects (see (usage)[#usage] below for examples) as well as a standalone CLI utility.

`CIGARCO` provides both a fine-grained classes to work it single alignments/queries and respective coordinate transformation requests, 
as well as for cases where there are multiple alignments and multiple coordinate transformation tasks.

#### Computational complexity
The main idea behind an efficient coordinate transformation comes from the concept of prefix sum arrays ([link](https://en.wikipedia.org/wiki/Prefix_sum)), 
which we compute (once on first access) for a query and target sequences given their CIGAR alignment, and which encode the number of query/target consumed positions at every distinct operation in the CIGAR string.
Such computation takes `O(n)` time, where `n` is the number of operations encoded in the CIGAR string.
The coordinate  transformation request is handled then via a look up for the supplied coordinate within a CIGAR string via a binary search in the respective prefix sum array.
This operation takes `O(log(n))` time. Then a matching number of consumed bases in the aligned sequenced is retrieved, and the overall coordinate transformation is computed (all steps beyond initial prefix sum array index lookup tak `O(1)` time).
The coordinate transformation requests are cached so the lookup of the exactly same coordinate takes `O(1)` time.

Overall, the complexity for a single query/target sequence alignment coordinate transformation takes `O(n)` on first invocation (computation of prefix sum arrays),
`O(log(n))` on second and subsequent invocation for coordinate values that have not been previously used, and `O(1)` for recurrent input values.   

#### Assumptions and ambiguous situations 
When handling multiple alignments at the same time, an assumption is made that every query sequence has exactly one alignment (i.e., no supplementary/secondary alignments are supported).
If multiple alignments for the same query sequence are provided -- the last one is used for computations.

Handling of coordinate transformation for reverse alignments, as well as for target -> query cases, are currently only supported at the CMapper API level, and not at the CLI application level.
From now on, when the target -> query transformation is considered all references to query and target shall be considered reversed, respectively.

If a requested source coordinate falls onto an insertion sequence in the query, the result of coordinate transformation will equal to the index of the latest matching position in the alignment prior to the insertion in question, or alignment start, if not such matching exists.
Coordinates for transformation have to fall in range between 0 and query length (as determined by the query consuming operations in the alignment CIGAR string ).
Coordinate `0` is always suitable for transformation, even with edge-case 0-length alignments (e.g., `0M`, while not real, is permitted by SAM specification).
When considering a reverse alignment and a query -> target transformation, `CIGARCO` assumes that the coordinate is given w.r.t. the original read, and not the aligned reverse complement.
When considering a reverse alignment and a target -> query transformation, coordinates are considered w.r.t. to regular target coordinate system.

#### Currently implemented functionality:
* transformation of query position into target coordinate system for both regular and reverse alignments 
* transformation of target position into query coordinate system for both regular and reverse alignments

For additional functionality that is coming in the future look at Future directions [section](#future-directions).  

### Installation:

    git clone https://github.com/aganezov/B2_CIGARCO.git
    cd B2_CIGARCO
    python setup.py install

### Usage

#### API
creating alignment object (immutable):
```python
alignment = Alignment(query_name="TR1", target_name="CHR1", start=3, cigar="8M7D6M2I2M11D7M", direction=True)
```

creating a mapper for an alignment object:
```python
mapper = CMapper(alignment)
```

transforming coordinates (query -> target):
```python
mapper.transform_coordinate(4)
7
mapper.transform_coordinate(13)
23
```

transforming coordinates (target -> query):
```python
mapper.transform_coordinate(4, direction='TQ')
1
mapper.transform_coordinate(40, direction='TQ')
21
```

#### CLI
```bash
>>> cigarco -a test/data/ex1_als.tsv -q test/data/ex1_qs.tsv
TR1 4   CHR1    7
TR2 0   CHR2    10
TR1 13  CHR1    23
TR2 10  CHR2    20
```

### Future direction


## Issues
Any issues shall be reported via an integrated GitHub [issue tracker](https://github.com/aganezov/B2_CIGARCO/issues).
For any complex questions please contact Sergey Aganezov who maintains the project at sergeyaganezojr[@]gmail.com.   

## Contributing

`CIGARCO` is developed with the TDD methodology, do for any future additions please consider creating a respective test suite to ensure that any and all added functionality is covered by tests.
TDD shall be implemented with already utilized test frameworks, which are `pytest` ([link](https://docs.pytest.org/en/stable/)) and `hypothesis` ([link](https://hypothesis.readthedocs.io/en/latest/)).


