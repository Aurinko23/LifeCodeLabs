import os
import pytest
from fasstq_filter import filter_fastq

def test_no_bound():
    try:
        filter_fastq("./input.fastq", "./no_bound_output.fastq")
        with open('./input.fastq', 'r') as f1, open('./no_bound_output.fastq', 'r') as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
        assert lines1 == lines2, "input file not same as output file"
    finally:
        if os.path.exists("./no_bound_output.fastq"):
            os.remove("./no_bound_output.fastq")


@pytest.mark.parametrize(
    "params",
    [
        ({"length_bounds": 1}),
        ({"gc_bounds": 1}),
        ({"quality_threshold": 101}),
    ]
)
def test_bound_all(params):
    try:
        filter_fastq("./input.fastq", "./no_bound_all.fastq", **params)
        with open('./no_bound_all.fastq', 'r') as f1:
            lines1 = f1.readlines()
        assert not lines1, "file not empty"
    finally:
        if os.path.exists("./no_bound_all.fastq"):
            os.remove("./no_bound_all.fastq")


@pytest.mark.parametrize(
    "params_int, params_tuple",
    [
        ({"length_bounds": 50}, {"length_bounds": (0, 50)}),
        ({"gc_bounds": 50}, {"gc_bounds": (0, 50)}),
    ]
)
def test_int_tuple(params_int, params_tuple):
    try:
        filter_fastq("./input.fastq", "./test_int.fastq", **params_int)
        filter_fastq("./input.fastq", "./test_tuple.fastq", **params_tuple)
        with open("./test_int.fastq", 'r') as f1, open("./test_tuple.fastq", 'r') as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
        assert lines1 == lines2, "int file not same as tuple file"
    finally:
        if os.path.exists("./test_int.fastq"):
            os.remove("./test_int.fastq")
        if os.path.exists("./test_tuple.fastq"):
            os.remove("./test_tuple.fastq")


def test_input_error():
    try:
        filter_fastq("./bad_input.fastq", "./bad_output_file.fastq")
        assert "never come here"
    except ValueError:
        pass

@pytest.mark.parametrize(
    "length_bounds, gc_bounds, quality_threshold",
    [
        (0, 50, -11),
    ]
)
def test_logger(length_bounds, gc_bounds, quality_threshold):
    try:
        filter_fastq("./input.fastq", "./outlog.fastq", length_bounds, gc_bounds, quality_threshold)
        with open('./error.log', 'r') as f1:
            lines1 = f1.readlines()
        assert lines1[-1] == "ERROR - quality_threshold=-11, can't be lower then 0, change to 0\n", "wrong log"
    finally:
        if os.path.exists("./outlog.fastq"):
            os.remove("./outlog.fastq")


def test_logger():
    try:
        filter_fastq("./input.fastq", "./outlog.fastq", 0, 50, -11)
        with open('./error.log', 'r') as f1:
            lines1 = f1.readlines()
        assert lines1[-1] == "ERROR - quality_threshold=-11, can't be lower then 0, change to 0\n", "wrong log"
    finally:
        if os.path.exists("./outlog.fastq"):
            os.remove("./outlog.fastq")