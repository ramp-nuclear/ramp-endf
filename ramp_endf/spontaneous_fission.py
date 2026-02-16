"""ENDF Spontaneous fission parsing module.

"""
from pathlib import Path
from typing import Sequence, Dict

from isotopes import ZAID
from .fission import parse_induced_fission, FissBranch, modlogger

SPFData = Dict[ZAID, Sequence[FissBranch]]


def parse_spontaneous_fission(*spf_paths: Path) -> SPFData:
    """Parse spontaneous fission data from the appropriate file"""

    modlogger.info("Adding spontaneous fission data")
    db = parse_induced_fission(*spf_paths)
    assert all([len(v) == 1 for v in db.values()])
    return {iso: val[0].branching for iso, val in db.items()}
