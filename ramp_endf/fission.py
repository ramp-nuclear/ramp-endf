"""Induced fission related tools to get data from ENDF.

"""
import itertools as it
import logging
from dataclasses import dataclass
from io import StringIO
from os import PathLike
from pathlib import Path
from typing import IO, Dict, Sequence, Tuple

import numpy as np
from endf.records import get_head_record, get_list_record
from isotopes import ZAID
from uncertainties import nominal_value, std_dev

from .evaluation import Evaluation, get_evaluations
from .util import ufloat

eV = float
FissBranch = Tuple[ZAID, float, float]

modlogger = logging.getLogger('endf.fission')


class FissionProductYields:
    """Independent and cumulative fission product yields.

    Parameters
    ----------
    ev_or_filename: Evaluation | PathLike
        ENDF fission product yield evaluation to read from. If given as a
        string, it is assumed to be the filename for the ENDF file.

    Attributes
    ----------
    cumulative : list of dict
        Cumulative yields for each tabulated energy. Each item in the list is a
        dictionary whose keys are nuclide names and values are cumulative
        yields. The i-th dictionary corresponds to the i-th incident neutron
        energy.
    energies : Iterable of float or None
        Energies at which fission product yields are tabulated.
    independent : list of dict
        Independent yields for each tabulated energy. Each item in the list is a
        dictionary whose keys are nuclide names and values are independent
        yields. The i-th dictionary corresponds to the i-th incident neutron
        energy.
    nuclide : dict
        Properties of the fissioning nuclide.

    Notes
    -----
    Neutron fission yields are typically not measured with a monoenergetic
    source of neutrons. As such, if the fission yields are given at, e.g.,
    0.0253 eV, one should interpret this as meaning that they are derived from a
    typical thermal reactor flux spectrum as opposed to a monoenergetic source
    at 0.0253 eV.

    """

    def __init__(self, ev_or_filename: Evaluation | PathLike):
        # Define function that can be used to read both independent and
        # cumulative yields
        def get_yields(file_obj: IO) -> tuple[np.ndarray, list[dict[ZAID, float]]]:
            # Determine number of energies
            n_energy = get_head_record(file_obj)[2]
            energies = np.zeros(n_energy)

            data = []
            for i in range(n_energy):
                # Determine i-th energy and number of products
                items, values = get_list_record(file_obj)
                energies[i] = items[0]
                n_products = items[5]

                # Get yields for i-th energy
                yields = {}
                for j in range(n_products):
                    Z, A = divmod(int(values[4 * j]), 1000)
                    isomeric_state = int(values[4 * j + 1])
                    zaid = ZAID(Z, A, isomeric_state)
                    yields[zaid] = ufloat(values[4 * j + 2], values[4 * j + 3])

                data.append(yields)

            return energies, data

        ev = ev_or_filename if isinstance(ev_or_filename, Evaluation) else Evaluation(ev_or_filename)

        # Assign basic nuclide properties
        self.zaid = ev.zaid
        self.nuclide = {
            'name': self.zaid.name,
            'atomic_number': ev.target['atomic_number'],
            'mass_number': ev.target['mass_number'],
            'isomeric_state': ev.target['isomeric_state']
        }

        # Read independent yields
        if (8, 454) in ev.section:
            file_obj = StringIO(ev.section[8, 454])
            self.energies, self.independent = get_yields(file_obj)

        # Read cumulative yields
        if (8, 459) in ev.section:
            file_obj = StringIO(ev.section[8, 459])
            energies, self.cumulative = get_yields(file_obj)
            assert np.all(energies == self.energies)

    @classmethod
    def from_endf(cls, ev_or_filename):
        """Generate fission product yield data from an ENDF evaluation
        
        Parameters
        ----------
        ev_or_filename : str or openmc.data.endf.Evaluation
            ENDF fission product yield evaluation to read from. If given as a
            string, it is assumed to be the filename for the ENDF file.
            
        Returns
        -------
        openmc.data.FissionProductYields
            Fission product yield data

        """
        return cls(ev_or_filename)


@dataclass(frozen=True, init=True)
class FissionYieldData:
    """Data object for fission yields

    energy - Inducing particle energy, in eV
    interpolation - Interpolation scheme. See ENDF-6 Formats Manual pg 22-27
    branching - Sequence of (target isotope, branching, branching error)

    """
    energy: eV
    interpolation: int
    branching: Sequence[FissBranch]


FissData = Dict[ZAID, Sequence[FissionYieldData]]


def parse_induced_fission(*nif_paths: Path, independent_yields = True) -> FissData:
    """parse induced fission data"""
    evals = it.chain.from_iterable(get_evaluations(path) for path in nif_paths)
    result = {}
    for evaluation in evals:
        data = FissionProductYields(evaluation)
        if independent_yields:
            fission_yield_data = data.independent
        else:
            fission_yield_data = data.cumulative
        result[data.zaid] = [
            FissionYieldData(ener, 0,
                             [(k, nominal_value(v), std_dev(v)) for k, v in yields.items()
                              if nominal_value(v) > 0.0])
            for ener, yields in zip(data.energies, fission_yield_data)]
    return result
