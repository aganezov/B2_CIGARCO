import logging
import os
from typing import List, Iterator
from unittest import mock

import hypothesis.strategies as st
import pytest
from hypothesis import given, settings

from cigarco.app import CigarcoApp, AlignmentsStreamer, QueryStreamer, create_cli_parser
from cigarco.mapping import CMapper, TransformationQuery, Alignment, TransformedResult
from test.test_cigar import decomposed_cigars
from test.test_manager import alignments as alignments_st


def test_alignment_streamer_setup_invalid_empty_constructor():
    with pytest.raises(TypeError):
        AlignmentsStreamer()


@given(iterable_data=st.lists(st.text()))
def test_alignment_streamer_setup(iterable_data):
    streamer = AlignmentsStreamer(iterable_data)
    assert streamer.alignments_str_source is iterable_data


@st.composite
def alignment_str_data_entries(draw):
    query_name = draw(st.text(min_size=1).filter(lambda e: "\t" not in e).filter(lambda e: len(e.strip()) > 0))
    target_name = draw(st.text(min_size=1).filter(lambda e: "\t" not in e).filter(lambda e: len(e.strip()) > 0))
    coordinate = str(draw(st.integers(min_value=0)))
    cigar_str = "".join(f"{cnt}{op}" for cnt, op in zip(*draw(decomposed_cigars())))
    return "\t".join([query_name, target_name, coordinate, cigar_str])


@given(alignment_str_data=st.lists(alignment_str_data_entries()))
def test_alignment_streamer_iteration(alignment_str_data):
    streamer = AlignmentsStreamer(alignment_str_data)
    alignments: Iterator[Alignment] = iter(streamer)
    result: List[Alignment] = [al for al in alignments]
    assert all(map(lambda al: isinstance(al, Alignment), result))
    assert len(result) == len(alignment_str_data)


@st.composite
def alignment_str_data_entries_too_few_entries(draw):
    alignment_str_entry = draw(alignment_str_data_entries())
    str_entries = alignment_str_entry.split("\t")
    removed_index = draw(st.integers(min_value=0, max_value=3))
    str_entries.pop(removed_index)
    return "\t".join(str_entries)


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       incorrect_alignment_str_data=alignment_str_data_entries_too_few_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()))
def test_alignment_streamer_iteration_invalid_entry_too_few_values_fail(alignment_str_data1, incorrect_alignment_str_data, alignment_str_data2):
    alignment_str_data = alignment_str_data1 + [incorrect_alignment_str_data] + alignment_str_data2
    streamer = AlignmentsStreamer(alignment_str_data)
    with pytest.raises(ValueError):
        [al for al in iter(streamer)]


@st.composite
def invalid_coordinate_str_entries(draw):
    invalid_numeric_value = str(draw(st.integers(max_value=-1)))
    non_numerical_value = draw(st.text().filter(lambda e: not e.isnumeric() and e != "-0"))
    selection = draw(st.integers(min_value=0, max_value=1))
    return [invalid_numeric_value, non_numerical_value][selection]


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       alignment_str_data=alignment_str_data_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()),
       invalid_coordinate_str_entry=invalid_coordinate_str_entries())
def test_alignment_streamer_iteration_invalid_entry_invalid_coordinate(alignment_str_data1, alignment_str_data, alignment_str_data2, invalid_coordinate_str_entry):
    alignment_str_data = alignment_str_data.split("\t")
    alignment_str_data[2] = invalid_coordinate_str_entry
    alignment_str_data = "\t".join(alignment_str_data)
    alignment_str_data = alignment_str_data1 + [alignment_str_data] + alignment_str_data2
    streamer = AlignmentsStreamer(alignment_str_data)
    with pytest.raises(ValueError):
        [al for al in iter(streamer)]


@st.composite
def invalid_cigar_strings(draw):
    decomposed_cigar = [(cnt, op) for cnt, op in zip(*draw(decomposed_cigars()))]
    if len(decomposed_cigar) == 0:
        return "M"
    changed_site = draw(st.integers(min_value=0, max_value=len(decomposed_cigar)))
    decomposed_cigar[changed_site] = (-1, "M")
    return "".join(f"{cnt}{op}" for cnt, op in decomposed_cigar)


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       alignment_str_data=alignment_str_data_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()),
       invalid_cigar=invalid_coordinate_str_entries())
def test_alignment_streamer_iteration_invalid_entry_invalid_cigar(alignment_str_data1, alignment_str_data, alignment_str_data2, invalid_cigar):
    alignment_str_data = alignment_str_data.split("\t")
    alignment_str_data[3] = invalid_cigar
    alignment_str_data = "\t".join(alignment_str_data)
    alignment_str_data = alignment_str_data1 + [alignment_str_data] + alignment_str_data2
    streamer = AlignmentsStreamer(alignment_str_data)
    with pytest.raises(ValueError):
        [al for al in iter(streamer)]


def test_query_streamer_setup_invalid_empty_constructor():
    with pytest.raises(TypeError):
        QueryStreamer()


@st.composite
def query_strings(draw, queries_names=None, max_coordinate=None):
    if queries_names is None:
        query_name = draw(st.text().filter(lambda e: "\t" not in e and len(e.strip()) > 0))
    else:
        query_name = st.one_of(queries_names)
    coordinate = str(draw(st.integers(min_value=0, max_value=1000000000 if max_coordinate is None else max_coordinate)))
    return "\t".join([query_name, coordinate])


@given(query_string_list=st.lists(query_strings()))
def test_query_streamer_setup(query_string_list):
    streamer = QueryStreamer(query_string_list)
    assert streamer.queries_str_source is query_string_list


@given(query_string_list=st.lists(query_strings()))
def test_query_streamer_iteration(query_string_list):
    streamer = QueryStreamer(query_string_list)
    queries = iter(streamer)
    result = [q for q in queries]
    assert len(result) == len(query_string_list)
    assert all(map(lambda e: isinstance(e, TransformationQuery), result))


@given(query_string_list1=st.lists(query_strings()),
       query_string=query_strings(),
       query_string_list2=st.lists(query_strings()))
def test_query_streamer_iteration_invalid_too_few_entries(query_string_list1, query_string, query_string_list2):
    query_string = query_string.split("\t")[0]
    query_strings = query_string_list1 + [query_string] + query_string_list2
    streamer = QueryStreamer(query_strings)
    with pytest.raises(ValueError):
        [q for q in streamer]


@given(query_string_list1=st.lists(query_strings()),
       query_string=query_strings(),
       query_string_list2=st.lists(query_strings()),
       invalid_coordinate_str=st.text().filter(lambda e: not e.isnumeric()))
def test_test_query_streamer_iteration_invalid_coordinate(query_string_list1, query_string, query_string_list2, invalid_coordinate_str):
    split_q_string = query_string.split("\t")
    query_string = split_q_string[0] + "\t" + invalid_coordinate_str
    query_strings = query_string_list1 + [query_string] + query_string_list2
    streamer = QueryStreamer(query_strings)
    with pytest.raises(ValueError):
        [q for q in streamer]


def test_cigarco_app_setup_invalid():
    with pytest.raises(TypeError):
        CigarcoApp()


@st.composite
def transformation_inputs(draw, max_al_size: int = 100000000, al_cnts: int = 1, max_q_cnt: int = 1):
    source_alignments = draw(st.lists(alignments_st(max_cigar_op_cnt=max_al_size), min_size=1, max_size=max(1, al_cnts)))
    # because of the fact that only one alignment is supported per query at this stage
    # we need to ensure that no two alignments share a query name
    unique_query_names = {al.query_name: al for al in source_alignments}
    source_alignments = list(unique_query_names.values())
    queries = []
    for alignment in source_alignments:
        query_length = CMapper(alignment).query_prefix_sums[-1]
        query_cnt = draw(st.integers(min_value=1, max_value=max(1, max_q_cnt)))
        for _ in range(query_cnt):
            coordinate = draw(st.integers(min_value=0, max_value=max(0, query_length - 1)))
            queries.append(TransformationQuery(alignment.query_name, coordinate))
    return source_alignments, queries


@given(transform_input=transformation_inputs())
def test_cigarco_app_setup(caplog, transform_input):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        source_alignments, transform_queries = transform_input
        app = CigarcoApp(source_alignments, transform_queries)
        assert app.fail_mode == 'R'
    assert len(caplog.text) > 0


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       alignment_str=alignment_str_data_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()))
def test_cigarco_app_setup_invalid_alignments_report(caplog, alignment_str_data1, alignment_str, alignment_str_data2):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        alignment_str_data = alignment_str.split("\t")
        alignment_str_data[2] = "-1"
        alignment_str = "\t".join(alignment_str_data)
        alignment_strs = alignment_str_data1 + [alignment_str] + alignment_str_data2
        CigarcoApp(AlignmentsStreamer(alignment_strs), [], fail_mode="R")
    assert caplog.text.count("ERROR") == 1


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       alignment_str=alignment_str_data_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()))
def test_cigarco_app_setup_invalid_alignments_ignore(caplog, alignment_str_data1, alignment_str, alignment_str_data2):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        alignment_str_data = alignment_str.split("\t")
        alignment_str_data[2] = "-1"
        alignment_str = "\t".join(alignment_str_data)
        alignment_strs = alignment_str_data1 + [alignment_str] + alignment_str_data2
        CigarcoApp(AlignmentsStreamer(alignment_strs), [], fail_mode="I")
    assert caplog.text.count("ERROR") == 0


@given(alignment_str_data1=st.lists(alignment_str_data_entries()),
       alignment_str=alignment_str_data_entries(),
       alignment_str_data2=st.lists(alignment_str_data_entries()))
def test_cigarco_app_setup_invalid_alignments_fail(caplog, alignment_str_data1, alignment_str, alignment_str_data2):
    with pytest.raises(SystemExit) as pytest_wrapper_r:
        caplog.clear()
        with caplog.at_level(logging.INFO):
            alignment_str_data = alignment_str.split("\t")
            alignment_str_data[2] = "-1"
            alignment_str = "\t".join(alignment_str_data)
            alignment_strs = alignment_str_data1 + [alignment_str] + alignment_str_data2
            CigarcoApp(AlignmentsStreamer(alignment_strs), [], fail_mode="F")
    assert caplog.text.count("CRITICAL") == 1
    assert pytest_wrapper_r.type == SystemExit
    assert pytest_wrapper_r.value.code == 1


@given(transformation_input=transformation_inputs(max_al_size=1000, al_cnts=100))
def test_cigarco_app_mapping_iteration(caplog, transformation_input):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        alignments, queries = transformation_input
        app = CigarcoApp(alignments, queries)
        results = [tr for tr in app.transformations_iter()]
        assert "Processing coordinate transformation queries" in caplog.text
        assert len(results) == len(queries)
        assert all(map(lambda tr: isinstance(tr[0], TransformationQuery), results))
        assert all(map(lambda tr: isinstance(tr[1], TransformedResult), results))
        assert all(map(lambda tr: tr[1].coordinate >= 0, results))


@given(transformation_input=transformation_inputs(max_al_size=1000, al_cnts=100))
def test_cigarco_app_mapping_iteration_invalid_query_setup_ignore(caplog, transformation_input):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        alignments, queries = transformation_input
        alignment: Alignment = alignments[0]
        q_name = alignment.query_name
        queries.append(TransformationQuery(q_name, coordinate=-1))
        app = CigarcoApp(alignments, queries, fail_mode="I")
        results = [tr for tr in app.transformations_iter()]
        assert "Processing coordinate transformation queries" in caplog.text
        assert len(results) == len(queries) - 1
        assert all(map(lambda tr: isinstance(tr, tuple), results))
        assert all(map(lambda tr: isinstance(tr[0], TransformationQuery), results))
        assert all(map(lambda tr: isinstance(tr[1], TransformedResult), results))
        assert all(map(lambda tr: tr[1].coordinate >= 0, results))


@given(transformation_input=transformation_inputs(max_al_size=1000, al_cnts=100))
def test_cigarco_app_mapping_iteration_invalid_query_setup_report(caplog, transformation_input):
    caplog.clear()
    with caplog.at_level(logging.INFO):
        alignments, queries = transformation_input
        alignment: Alignment = alignments[0]
        q_name = alignment.query_name
        queries.append(TransformationQuery(q_name, coordinate=-1))
        app = CigarcoApp(alignments, queries, fail_mode="R")
        results = [tr for tr in app.transformations_iter()]
        assert caplog.text.count("ERROR") == 1
        assert len(results) == len(queries) - 1
        assert all(map(lambda tr: isinstance(tr, tuple), results))
        assert all(map(lambda tr: isinstance(tr[0], TransformationQuery), results))
        assert all(map(lambda tr: isinstance(tr[1], TransformedResult), results))
        assert all(map(lambda tr: tr[1].coordinate >= 0, results))


@settings(deadline=None)
@given(transformation_input=transformation_inputs(max_al_size=1000, al_cnts=100))
def test_cigarco_app_mapping_iteration_invalid_query_setup_fail(caplog, transformation_input):
    with pytest.raises(SystemExit) as pytest_wrapper_r:
        caplog.clear()
        with caplog.at_level(logging.INFO):
            alignments, queries = transformation_input
            alignment: Alignment = alignments[0]
            q_name = alignment.query_name
            queries.append(TransformationQuery(q_name, coordinate=-1))
            app = CigarcoApp(alignments, queries, fail_mode="F")
            [tr for tr in app.transformations_iter()]
    assert caplog.text.count("CRITICAL") == 1
    assert pytest_wrapper_r.type == SystemExit
    assert pytest_wrapper_r.value.code == 1


def test_cigarco_app_default_parser_setup():
    parser = create_cli_parser()
    args = parser.parse_args(["-a" "test/data/ex1_als.tsv", "-q" "test/data/ex1_qs.tsv"])
    assert args.alignments
    assert args.queries
    assert args.error_mode == "R"
    assert args.log_level == 20
    assert args.output


def test_cigarco_app_execution():
    from cigarco import app
    with mock.patch.object(app, "__name__", "__main__"):
        out_file_path = "test/data/ex1_out.tsv"
        if os.path.exists(out_file_path):
            os.remove(out_file_path)
        import sys
        sys.argv = ["python",
                    "-a" "test/data/ex1_als.tsv",
                    "-q" "test/data/ex1_qs.tsv",
                    "-o", out_file_path]
        app.execute_script()
        with open(out_file_path, "rt") as source:
            produced_data = [l.strip() for l in source.readlines()]
            assert len(produced_data) == 4
            assert "\t".join(["TR1", "4", "CHR1", "7"]) in produced_data
            assert "\t".join(["TR2", "0", "CHR2", "10"]) in produced_data
            assert "\t".join(["TR1", "13", "CHR1", "23"]) in produced_data
            assert "\t".join(["TR2", "10", "CHR2", "20"]) in produced_data
        if os.path.exists(out_file_path):
            os.remove(out_file_path)
