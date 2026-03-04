"""Utility functions for reading ENDF6 formats"""

import logging
import math
import re
from typing import Any, Generator, Iterable, Sequence, TextIO, Union

import numpy as np
from isotopes import ZAID
from uncertainties import UFloat
from uncertainties import ufloat as raw_ufloat

from .units import CMBarnPerSecond, Second


def _za_to_isotope(za: float, m: int) -> ZAID:
    z, a = divmod(np.round(za), 1000)
    return ZAID(z, a, m)


_FORTRAN_PAT = re.compile(r"(?P<mantisa>\s*[+-]?[\d.]+)(?P<exp>[+-]\d{1,2})")
modlogger = logging.getLogger("endf.util")


def ufloat(mu: float, sigma: float) -> UFloat | float:
    """Returns a ufloat if the standard deviation isn't 0, otherwise returns the float"""
    return raw_ufloat(mu, sigma) if sigma != 0.0 else mu


def cmplx_float(s: str):
    """Parse a float that can be written in obnoxious fortran ways.

    Paremeters
    ----------
    s - string to convert

    Returns
    -------

    float version of that stupidly formatted string

    """

    try:
        return float(s)
    except ValueError as e:
        m = re.match(_FORTRAN_PAT, s)
        if not m:
            raise e
        return float(m.group("mantisa")) * (10 ** int(m.group("exp")))


def _split_line(header: str) -> Iterable[str]:
    return (header[i:j] for i, j in zip(range(0, 66, 11), range(11, 77, 11)))


def _lineparse(_lines: Sequence[str]) -> Generator[float, Any, Any]:
    for line in _lines:
        for val in _split_line(line):
            yield cmplx_float(val)


# noinspection PyPep8Naming
def parse_section(of: TextIO, MF: int, MT: int, nuclide: int = 0):
    """Parse an MT section from a file.

    Parameters
    ----------
    of = Open file IO
    MF = File number to read. Usually 8 for decay data in ENDF6 format.
    MT = Card number to read. Usually 8 for decay data in ENDF6 format.
    nuclide = Nuclide identifier in the library. If not given, assume there is
              only one within the file IO.

    Returns
    -------

    """

    pat = re.compile(
        "^(" + 66 * "." + (f")(?:{nuclide:4d})" if nuclide else ")(?:....)") + f"(?:{MF:2d}{MT:3d})", re.MULTILINE
    )
    # noinspection PyTypeChecker
    return re.findall(pat, of.read())


def parse_list(lines: Sequence[str], num_vector: int, num_total: int):
    """Parse a list record from the ENDF section.

    Parameters
    ----------
    lines - ENDF section lines, starts at where the list starts.
    num_of_values - Number of floating values to read.

    Returns
    -------

    A tuple of vectors and the remaining lines after reading out the
    relevant list
    """

    modlogger.debug("Reading a list for %d vectors, total of %d numbers", num_vector, num_total)
    data = np.zeros((num_total // num_vector, num_vector), dtype="float32")
    nlines = int(np.ceil(num_total / 6))
    modlogger.debug("Reading %d lines to get data...", nlines)
    rem = lines[nlines:]
    for n, v in zip(range(num_total), _lineparse(lines[:nlines])):
        data.flat[n] = v
    return tuple([data[:, i] for i in range(num_vector)] + [rem])


def parse_header(header: str) -> Generator[Union[float, int], Any, Any]:
    """Parses a CONT or HEAD record.

    These are fortran typed by:
    2E11.0, 4I11. The actual full record includes the MAT-MF-MT-NS tail, but
    since these are cleaned by the parse_section function, this function
    assumes they are not there.

    TODO: Make this safe when numbers touch each other.

    Parameters
    ----------
    header - The record to parse

    Returns
    -------

    The parsed data, first two floats and then 4 integers.
    """
    return (cmplx_float(v) if i < 2 else int(v) for i, v in enumerate(_split_line(header)))


def halflife_to_rate(halflife: Second, branching: float = 1.0) -> CMBarnPerSecond:
    """Get the decay rate of an isotope from its halflife.

    Parameters
    ----------
    halflife - Halflife of this mode of decay
    branching - Branching factor to this mode of decay

    Returns
    -------
    The rate at which the isotope decays in this mode.

    """
    return math.log(2) * branching / halflife
