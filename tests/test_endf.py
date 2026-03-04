#!/usr/bin/env python3
"""Theses tests should make sure that endf reads properly from sample files."""

import json
from math import isclose
from os.path import join
from pathlib import Path

import pytest
from networkx.readwrite import json_graph
from ramp_core.json import RampJSONEncoder

try:
    from batman.graphs.decay import DecayGraph
except ImportError:
    DecayGraph = None

from isotopes import (
    I135,
    U238,
    Cs135,
    He4,
    Th234,
    Xe135,
    Xe135m,
)

from ramp_endf import parse_decay_processes, parse_spontaneous_fission
from ramp_endf.modes import ALPHA, SPF

test_data_dir = Path(__file__).parent / "test_data"
slow = pytest.mark.slow
needs_batman = pytest.mark.skipif(reason="Can't import batman", condition=DecayGraph is None)


def needs_data(*fnames):
    paths = [test_data_dir / fname for fname in fnames]

    def _data_taker(f):
        not_located = [path for path in paths if not path.exists()]
        return pytest.mark.skipif(
            reason=f"Can't locate {[str(v) for v in not_located]}. Use the download script in {str(test_data_dir)}",
            condition=not_located,
        )(f)

    return _data_taker


# noinspection PyPep8Naming
@needs_data("dec-053_I_135.endf")
def test_I135_file_read_has_certain_known_properties():
    f = Path(join(test_data_dir, "dec-053_I_135.endf"))
    processes = parse_decay_processes(f)
    iso = processes[0].parent
    assert iso == I135
    assert isclose(processes[0].halflife, 2.365200e4)
    assert processes[0].targets == {Xe135}
    assert processes[1].targets == {Xe135m}


# noinspection PyPep8Naming
@needs_data("dec-054_Xe_135.endf")
def test_Xe135_file_read_has_certain_known_properties():
    f = Path(join(test_data_dir, "dec-054_Xe_135.endf"))
    (process,) = parse_decay_processes(f)
    assert process.parent == Xe135
    assert process.targets == {Cs135}


@needs_data("sfy-092_U_238.endf", "dec-092_U_238.endf")
@needs_batman
def test_U238_file_read_has_certain_known_properties():
    spf_file = Path(join(test_data_dir, "sfy-092_U_238.endf"))
    decay_file = Path(join(test_data_dir, "dec-092_U_238.endf"))
    db = parse_spontaneous_fission(spf_file)
    assert db.keys() == {U238}
    processes = parse_decay_processes(decay_file, spf_db=db)
    iso = processes[0].parent
    assert iso == U238
    for process in processes:
        assert isclose(process.halflife, 1.40999e17)
        assert process.mode[0] in {ALPHA, SPF}
        if process.mode[0] == ALPHA:
            assert process.targets == {He4, Th234}
            assert isclose(process.energy, 4.2697e6)
            assert isclose(process.energy_err, 2.9e3)
            assert isclose(process.fraction, 1.0)
            assert isclose(process.fraction_err, 0.0)
        elif process.mode[0] == SPF:
            assert isclose(process.energy, 1.736e8)
            assert isclose(process.energy_err, 3.4e6)


@pytest.mark.regression
@needs_data("sfy-092_U_238.endf", "dec-092_U_238.endf")
@needs_batman
def test_U238_file_read_graph_by_regression(data_regression):
    spf_file = Path(join(test_data_dir, "sfy-092_U_238.endf"))
    decay_file = Path(join(test_data_dir, "dec-092_U_238.endf"))
    db = parse_spontaneous_fission(spf_file)
    processes = parse_decay_processes(decay_file, spf_db=db)
    iso = processes[0].parent
    g = DecayGraph()
    for process in processes:
        g.add_edge_from_process(process)
    data = json_graph.node_link_data(g)
    data_regression.check(json.dumps(data, cls=RampJSONEncoder))


if __name__ == "__main__":
    frmt = "%(filename)s - %(lineno)d - %(asctime)s - %(levelname)s: %(message)s"
    import logging
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Console utility to ease my life")
    parser.add_argument("f", help="File to read from")
    parser.add_argument("spf", default=None, help="Spontaneous fission file")
    parser.add_argument("--log", default="test.log", help="Logging file")
    parser.add_argument("--logmode", default="w", help="Logging mode")
    parser.add_argument("--format", default=frmt, help="Logging format")
    parser.add_argument("-d", "--debug", action="store_true", help="Flag to enable log debug mode")

    args = parser.parse_args()

    level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(filename=args.log, filemode=args.logmode, format=args.format, level=level)
    spf_data = {}
    if args.spf:
        logging.info("Starting up specific spontaneous fission parsing for %s", args.spf)
        spf_data.update(parse_spontaneous_fission(args.spf))
    logging.info("Starting up specific file parsing for %s...", args.f)
    processes = parse_decay_processes(args.f, spf_db=spf_data)
    par = processes[0].parent
    logging.info("Read information for Isotope %s", par)
    for process in processes:
        logging.info(f"Found a decay branch that leads to {process.targets}")
    else:
        logging.info("Succeeded.")
    logging.info("Printing log file to screen...")
    with open(args.log, "r") as fo:
        print(fo.read())
