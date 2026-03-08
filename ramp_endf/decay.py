#!/usr/bin/env python3
"""Module for reading decay data from ENDF files.

This uses the dec files, that contain decay data for card MT=451

"""

import itertools as it
from collections import Counter
from dataclasses import dataclass, field
from io import StringIO
from math import log
from pathlib import Path
from typing import IO, Any, Dict, FrozenSet, List, Optional, Sequence, Tuple, Type, TypeVar, Union

from endf.records import get_head_record, get_list_record, get_tab1_record
from isotopes import H1, ZAID, He4
from multipledispatch import dispatch
from uncertainties import UFloat, nominal_value, std_dev

from .evaluation import Evaluation, get_evaluations
from .modes import ALPHA, BETA_M, BETA_P, ISOMERIC_TRANSITION, NEUTRON, PROTON, SPF
from .rad_types import SF, e_m, n, p, x, α, β_m, β_p, γ
from .units import EV, CMBarnPerSecond, Second
from .util import halflife_to_rate, ufloat

__all__ = ["parse_decay_processes", "DecayProcess", "Decay"]


BOTH = "both"
CONTINUOUS = "continuous"
DISCRETE = "discrete"

eV = float

Sink = ZAID(0, 0, 0)

_RTYP_DECAY_MODES = {
    1: BETA_M,
    2: BETA_P,
    3: ISOMERIC_TRANSITION,
    4: ALPHA,
    5: NEUTRON,
    6: SPF,
    7: PROTON,
}
# Gives name and (change in Z, change in A) resulting from decay
_DECAY_MODES_CHANGES = {
    BETA_M: (1, 0),
    BETA_P: (-1, 0),
    ISOMERIC_TRANSITION: (0, 0),
    ALPHA: (-2, -4),
    NEUTRON: (0, -1),
    PROTON: (-1, -1),
}

_DECAY_MODES_PRODUCTS = {
    PROTON: H1,
    ALPHA: He4,
}
_RADIATION_TYPES = {
    0: γ,
    1: β_m,
    2: β_p,
    4: α,
    5: n,
    6: SF,
    7: p,
    8: e_m,
    9: x,
}
T = TypeVar("T", bound="Decay")


class TransitionType:
    Allowed = 0
    FirstForbidden = 2
    SecondForbidden = 3


SPFData = Dict[ZAID, Sequence[Tuple[ZAID, float, float]]]


@dataclass(init=True, frozen=True, unsafe_hash=True)
class DecayProcess:
    """Class meant to define the data in a radioactive decay process."""

    parent: ZAID
    _target: Union[Tuple[Tuple[ZAID, float], ...], ZAID]
    halflife: Second
    mode: Sequence[str]
    fraction: float = 1.0
    fraction_err: float = 0.0
    energy: EV = 0.0
    energy_err: EV = 0.0
    spectra: Dict[Any, Any] = field(default_factory=dict)

    @property
    def target_branching(self) -> Dict[ZAID, float]:
        """Isotope branching mapping"""
        if isinstance(self._target, ZAID):
            return {self._target: 1.0}
        return dict(self._target)

    @property
    def targets(self) -> FrozenSet[ZAID]:
        """Isotopes the decay can create"""
        return frozenset(self.target_branching.keys())

    @property
    def rate(self) -> CMBarnPerSecond:
        """Returns reaction rate per parent nucleus."""
        return halflife_to_rate(self.halflife, branching=self.fraction)

    @classmethod
    def from_other(cls, other: "DecayProcess", **kwargs) -> "DecayProcess":
        """Create a copy of a reaction with changes.

        Used mostly for debugging, to redefine processes where necessary.

        Parameters
        ----------
        other - A different reaction

        Returns
        -------

        Decay process that is mostly similar but with some changes.

        """

        for key in [
            "parent",
            "_target",
            "halflife",
            "mode",
            "fraction",
            "fraction_err",
            "energy",
            "energy_err",
            "spectra",
        ]:
            kwargs[key] = kwargs.get(key, other.__getattribute__(key))
        return cls(**kwargs)


def get_decay_modes(value):
    """Return sequence of decay modes given an ENDF RTYP value.

    Parameters
    ----------
    value : float
        ENDF definition of sequence of decay modes

    Returns
    -------
    Sequence[str]
        Successive decays, e.g. ('beta-', 'neutron')

    """
    return tuple(_RTYP_DECAY_MODES[int(v)] for v in str(value).strip("0").replace(".", ""))


@dataclass
class DecayMode:
    """Radioactive decay mode.

    Parameters
    ----------
    parent : ZAID
        Parent decaying nuclide
    modes : Sequence[str]
        Successive decay modes
    daughter_state : int
        Metastable state of the daughter nuclide
    energy : UFloat
        Total decay energy in eV available in the decay process.
    branching_ratio : UFloat
        Fraction of the decay of the parent nuclide which proceeds by this mode.
    """

    parent: ZAID
    modes: Sequence[str]
    daughter_state: int
    energy: UFloat
    branching_ratio: UFloat

    def _increment_iso(self, iso: ZAID, z: int, a: int) -> ZAID:
        return ZAID(iso.Z + z, iso.A + a, self.daughter_state)

    def target_branching(self, spf_db: Optional[SPFData] = None) -> Tuple[Tuple[ZAID, float], ...]:
        """Branching as an immutable sequence of pairs"""
        counter = Counter()
        iso = self.parent
        for mode in self.modes:
            if product := _DECAY_MODES_PRODUCTS.get(mode):
                counter.update([product])
            if changes := _DECAY_MODES_CHANGES.get(mode):
                iso = self._increment_iso(iso, *changes)
            if mode == SPF:
                try:
                    counter.update({target: yield_ for target, yield_, _ in spf_db[iso]})
                except KeyError:
                    Warning(f"{iso} not in spf_db: {spf_db.keys()}")
                iso = None
        if iso:
            counter.update([iso])
        return tuple(counter.items())

    @property
    def daughter(self) -> ZAID:
        """Returns the daughter ZAID of this mode"""
        iso = self.parent
        for mode in self.modes:
            if changes := _DECAY_MODES_CHANGES.get(mode):
                iso = self._increment_iso(iso, *changes, self.daughter_state)
            else:
                return Sink
        return iso

    def __repr__(self) -> str:
        return "<DecayMode: ({}), {} -> {}, {}>".format(
            ",".join(self.modes), self.parent, self.daughter, self.branching_ratio
        )


def _floatify(v) -> float:
    return v if type(v) is float else v.item()


class Decay:
    """Radioactive decay data.

    Parameters
    ----------
    ev_or_filename : str of openmc.data.endf.Evaluation
        ENDF radioactive decay data evaluation to read from. If given as a
        string, it is assumed to be the filename for the ENDF file.

    Attributes
    ----------

    average_energies : dict
        Average decay energies in eV of each type of radiation for decay heat
        applications.
    decay_constant : UFloat
        Decay constant in inverse seconds.
    half_life : UFloat
        Half-life of the decay in seconds.
    modes : list
        Decay mode information for each mode of decay.
    nuclide : dict
        Dictionary describing decaying nuclide with keys 'name',
        'excited_state', 'mass', 'stable', 'spin', and 'parity'.
    spectra : dict
        Resulting radiation spectra for each radiation type.

    """

    def __init__(self, ev_or_filename: Evaluation | Path | IO):
        # Get evaluation if str is passed
        if isinstance(ev_or_filename, Evaluation):
            ev = ev_or_filename
        else:
            ev = Evaluation(ev_or_filename)

        file_obj = StringIO(ev.section[8, 457])

        self.nuclide = {}
        self.modes = []
        self.spectra = {}
        self.average_energies = {}

        # Get head record
        items = get_head_record(file_obj)
        z, a = divmod(items[0], 1000)
        metastable = items[3]
        self.nuclide["atomic_number"] = z
        self.nuclide["mass_number"] = a
        self.nuclide["isomeric_state"] = metastable
        self.zaid = ZAID(z, a, metastable)
        self.nuclide["name"] = self.zaid.name
        self.nuclide["mass"] = items[1]  # AWR
        self.nuclide["excited_state"] = items[2]  # State of the original nuclide
        self.nuclide["stable"] = items[4] == 1  # Nucleus stability flag

        # Determine if radioactive/stable
        if not self.nuclide["stable"]:
            nsp = items[5]  # Number of radiation types

            # Half-life and decay energies
            items, values = get_list_record(file_obj)
            self.half_life = ufloat(_floatify(items[0]), _floatify(items[1]))
            nc = items[4] // 2
            pairs = [v for v in zip(values[::2], values[1::2])]
            ex = self.average_energies
            ex["light"] = ufloat(*pairs[0])
            ex["electromagnetic"] = ufloat(*pairs[1])
            ex["heavy"] = ufloat(*pairs[2])
            if nc == 17:
                ex["beta-"] = ufloat(*pairs[3])
                ex["beta+"] = ufloat(*pairs[4])
                ex["auger"] = ufloat(*pairs[5])
                ex["conversion"] = ufloat(*pairs[6])
                ex["gamma"] = ufloat(*pairs[7])
                ex["xray"] = ufloat(*pairs[8])
                ex["Bremsstrahlung"] = ufloat(*pairs[9])
                ex["annihilation"] = ufloat(*pairs[10])
                ex["alpha"] = ufloat(*pairs[11])
                ex["recoil"] = ufloat(*pairs[12])
                ex["SF"] = ufloat(*pairs[13])
                ex["neutron"] = ufloat(*pairs[14])
                ex["proton"] = ufloat(*pairs[15])
                ex["neutrino"] = ufloat(*pairs[16])

            items, values = get_list_record(file_obj)
            spin = items[0]
            if spin == -77.777:
                self.nuclide["spin"] = None
            else:
                self.nuclide["spin"] = spin
            self.nuclide["parity"] = items[1]  # Parity of the nuclide

            # Decay mode information
            n_modes = items[5]  # Number of decay modes
            for i in range(n_modes):
                decay_type = get_decay_modes(values[6 * i])
                isomeric_state = int(values[6 * i + 1])
                energy = ufloat(*map(_floatify, values[6 * i + 2 : 6 * i + 4]))
                branching_ratio = ufloat(*map(_floatify, values[6 * i + 4 : 6 * (i + 1)]))

                mode = DecayMode(self.zaid, decay_type, isomeric_state, energy, branching_ratio)
                self.modes.append(mode)

            discrete_type = {
                0.0: None,
                1.0: TransitionType.Allowed,
                2.0: TransitionType.FirstForbidden,
                3.0: TransitionType.SecondForbidden,
            }

            # Read spectra
            for i in range(nsp):
                spectrum = {}

                items, values = get_list_record(file_obj)
                # Decay radiation type
                spectrum["type"] = _RADIATION_TYPES[items[1]]
                # Continuous spectrum flag
                spectrum["continuous_flag"] = {0: "discrete", 1: "continuous", 2: "both"}[items[2]]
                spectrum["discrete_normalization"] = ufloat(*values[0:2])
                spectrum["energy_average"] = ufloat(*values[2:4])
                spectrum["continuous_normalization"] = ufloat(*values[4:6])

                ner = items[5]  # Number of tabulated discrete energies

                if not spectrum["continuous_flag"] == "continuous":
                    # Information about discrete spectrum
                    spectrum["discrete"] = []
                    for j in range(ner):
                        items, values = get_list_record(file_obj)
                        di = {
                            "energy": ufloat(*items[0:2]),
                            "from_mode": get_decay_modes(values[0]),
                            "type": discrete_type[values[1]],
                            "intensity": ufloat(*values[2:4]),
                        }
                        if spectrum["type"] == β_p:
                            di["positron_intensity"] = ufloat(*values[4:6])
                        elif spectrum["type"] == γ:
                            if len(values) >= 6:
                                di["internal_pair"] = ufloat(*values[4:6])
                            if len(values) >= 8:
                                di["total_internal_conversion"] = ufloat(*values[6:8])
                            if len(values) == 12:
                                di["k_shell_conversion"] = ufloat(*values[8:10])
                                di["l_shell_conversion"] = ufloat(*values[10:12])
                        spectrum["discrete"].append(di)

                if not spectrum["continuous_flag"] == "discrete":
                    # Read continuous spectrum
                    ci = {}
                    params, ci["probability"] = get_tab1_record(file_obj)
                    ci["type"] = get_decay_modes(params[0])

                    # Read covariance (Ek, Fk) table
                    lcov = params[3]
                    if lcov != 0:
                        items, values = get_list_record(file_obj)
                        ci["covariance_lb"] = items[3]
                        ci["covariance"] = zip(values[0::2], values[1::2])

                    spectrum["continuous"] = ci

                # Add spectrum to dictionary
                self.spectra[spectrum["type"]] = spectrum

        else:
            _ = get_list_record(file_obj)
            items, values = get_list_record(file_obj)
            self.nuclide["spin"] = items[0]
            self.nuclide["parity"] = items[1]
            self.half_life = ufloat(float("inf"), float("inf"))

    @property
    def decay_constant(self) -> float:
        """Natural exponent decay rate equivalent of half life"""
        return log(2.0) / self.half_life

    @classmethod
    def from_endf(cls: Type[T], ev_or_filename: Evaluation | Path | IO) -> T:
        """Generate radioactive decay data from an ENDF evaluation

        Parameters
        ----------
        spf_db: SPFData | None
            Spontaneous fission additional data. Defaults to no data.
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF radioactive decay data evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.

        Returns
        -------
        openmc.data.Decay
            Radioactive decay data

        """
        return cls(ev_or_filename)

    def decay_processes(self, spf_db: Optional[SPFData] = None) -> List[DecayProcess]:
        spf_db = spf_db or {}
        return [
            DecayProcess(
                mode.parent,
                mode.target_branching(spf_db=spf_db),
                nominal_value(self.half_life),
                mode.modes,
                fraction=nominal_value(mode.branching_ratio),
                fraction_err=std_dev(mode.branching_ratio),
                energy=nominal_value(mode.energy),
                energy_err=std_dev(mode.energy),
            )
            for mode in self.modes
        ]


@dispatch(Decay)
def _cast_decays(dec: Decay) -> tuple[Decay]:
    return (dec,)


@dispatch(Evaluation)
def _cast_decays(ev: Evaluation) -> tuple[Decay]:
    return (Decay.from_endf(ev),)


@dispatch(object)
def _cast_decays(path: Path) -> tuple[Decay, ...]:
    return tuple(Decay.from_endf(ev) for ev in get_evaluations(path))


def parse_decays(*file_or_ev: Path | Evaluation | Decay) -> Sequence[Decay]:
    """parse Decay from decay files"""
    return list(it.chain.from_iterable(_cast_decays(d) for d in file_or_ev))


def parse_decay_processes(*file_or_ev: Path | Evaluation | Decay, spf_db=Optional[SPFData]) -> Sequence[DecayProcess]:
    """parse decay processes from decay files"""
    decays = parse_decays(*file_or_ev)
    return list(it.chain.from_iterable(decay.decay_processes(spf_db=spf_db) for decay in decays))
