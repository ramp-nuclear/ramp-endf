from io import StringIO
from os import PathLike
from pathlib import PurePath
from typing import IO, Iterable

from endf.records import get_cont_record, get_head_record, get_text_record, int_endf
from isotopes import ZAID

_LIBRARY = {0: 'ENDF/B', 1: 'ENDF/A', 2: 'JEFF', 3: 'EFF',
            4: 'ENDF/B High Energy', 5: 'CENDL', 6: 'JENDL',
            17: 'TENDL', 18: 'ROSFOND', 21: 'SG-21', 31: 'INDL/V',
            32: 'INDL/A', 33: 'FENDL', 34: 'IRDF', 35: 'BROND',
            36: 'INGDB-90', 37: 'FENDL/A', 41: 'BROND'}

_SUBLIBRARY = {
    0: 'Photo-nuclear data',
    1: 'Photo-induced fission product yields',
    3: 'Photo-atomic data',
    4: 'Radioactive decay data',
    5: 'Spontaneous fission product yields',
    6: 'Atomic relaxation data',
    10: 'Incident-neutron data',
    11: 'Neutron-induced fission product yields',
    12: 'Thermal neutron scattering data',
    19: 'Neutron standards',
    113: 'Electro-atomic data',
    10010: 'Incident-proton data',
    10011: 'Proton-induced fission product yields',
    10020: 'Incident-deuteron data',
    10030: 'Incident-triton data',
    20030: 'Incident-helion (3He) data',
    20040: 'Incident-alpha data'
}


class Evaluation:
    """ENDF material evaluation with multiple files/sections

    Parameters
    ----------
    filename_or_obj : str or file-like
        Path to ENDF file to read or an open file positioned at the start of an
        ENDF material

    Attributes
    ----------
    info : dict
        Miscellaneous information about the evaluation.
    target : dict
        Information about the target material, such as its mass, isomeric state,
        whether it's stable, and whether it's fissionable.
    projectile : dict
        Information about the projectile such as its mass.
    reaction_list : list of 4-tuples
        List of sections in the evaluation. The entries of the tuples are the
        file (MF), section (MT), number of records (NC), and modification
        indicator (MOD).

    """

    def __init__(self, filename_or_obj: PathLike | IO):
        if isinstance(filename_or_obj, (str, PurePath)):
            fh = open(str(filename_or_obj), 'r')
        else:
            fh = filename_or_obj
        self.section = {}
        self.info = {}
        self.target = {}
        self.projectile = {}
        self.reaction_list = []

        mf = 0
        # Determine MAT number for this evaluation
        while mf == 0:
            position = fh.tell()
            line = fh.readline()
            mf = int_endf(line[70:72])
        # These are actually safe because the while loop is always entered at least once.
        self.material = int_endf(line[66:70])
        fh.seek(position)

        while True:
            # Find next section
            while True:
                position = fh.tell()
                line = fh.readline()
                mat = int_endf(line[66:70])
                mf = int_endf(line[70:72])
                mt = int_endf(line[72:75])
                if mt > 0 or mat == 0:
                    fh.seek(position)
                    break

            # If end of material reached, exit loop
            if mat == 0:
                fh.readline()
                break

            section_data = ''
            while True:
                line = fh.readline()
                if line[72:75] == '  0':
                    break
                else:
                    section_data += line
            self.section[mf, mt] = section_data

        self._read_header()

    def __repr__(self):
        if 'zsymam' in self.target:
            name = self.target['zsymam'].replace(' ', '')
        else:
            name = 'Unknown'
        return '<{} for {} {}>'.format(self.info['sublibrary'], name,
                                       self.info['library'])

    def _read_header(self):
        file_obj = StringIO(self.section[1, 451])

        # Information about target/projectile
        items = get_head_record(file_obj)
        Z, A = divmod(items[0], 1000)
        self.target['atomic_number'] = Z
        self.target['mass_number'] = A
        self.target['mass'] = items[1]
        self._LRP = items[2]
        self.target['fissionable'] = (items[3] == 1)
        try:
            library = _LIBRARY[items[4]]
        except KeyError:
            library = 'Unknown'
        self.info['modification'] = items[5]

        # Control record 1
        items = get_cont_record(file_obj)
        self.target['excitation_energy'] = items[0]
        self.target['stable'] = (int(items[1]) == 0)
        self.target['state'] = items[2]
        self.target['isomeric_state'] = m = items[3]
        self.zaid = ZAID(Z, A, m)
        self.info['format'] = items[5]
        assert self.info['format'] == 6

        # Set correct excited state for Am242_m1, which is wrong in ENDF/B-VII.1
        if Z == 95 and A == 242 and m == 1:
            self.target['state'] = 2

        # Control record 2
        items = get_cont_record(file_obj)
        self.projectile['mass'] = items[0]
        self.info['energy_max'] = items[1]
        library_release = items[2]
        self.info['sublibrary'] = _SUBLIBRARY[items[4]]
        library_version = items[5]
        self.info['library'] = (library, library_version, library_release)

        # Control record 3
        items = get_cont_record(file_obj)
        self.target['temperature'] = items[0]
        self.info['derived'] = (items[2] > 0)
        NWD = items[4]
        NXC = items[5]

        # Text records
        text = [get_text_record(file_obj) for i in range(NWD)]
        if len(text) >= 5:
            self.target['zsymam'] = text[0][0:11]
            self.info['laboratory'] = text[0][11:22]
            self.info['date'] = text[0][22:32]
            self.info['author'] = text[0][32:66]
            self.info['reference'] = text[1][1:22]
            self.info['date_distribution'] = text[1][22:32]
            self.info['date_release'] = text[1][33:43]
            self.info['date_entry'] = text[1][55:63]
            self.info['identifier'] = text[2:5]
            self.info['description'] = text[5:]

        # File numbers, reaction designations, and number of records
        for i in range(NXC):
            _, _, mf, mt, nc, mod = get_cont_record(file_obj, skip_c=True)
            self.reaction_list.append((mf, mt, nc, mod))


def get_evaluations(filename: PathLike) -> Iterable[Evaluation]:
    """Return a list of all evaluations within an ENDF file.

    Parameters
    ----------
    filename : str
        Path to ENDF-6 formatted file

    Returns
    -------
    list
        A list of :class:`openmc.data.endf.Evaluation` instances.

    """
    with open(str(filename), 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if not line or line[66:70] == '  -1':
                break
            fh.seek(pos)
            yield Evaluation(fh)
