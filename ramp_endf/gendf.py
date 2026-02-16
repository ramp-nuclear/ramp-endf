from os import PathLike
from pathlib import PurePath
from io import StringIO
from typing import IO

from endf.records import get_list_record, get_cont_record
import itertools as it


class GENDF:
    def __init__(self, filename_or_obj: PathLike | IO):
        if isinstance(filename_or_obj, (str, PurePath)):
            fh = open(str(filename_or_obj), 'r')
        else:
            fh = filename_or_obj
        self.section = {}
        self.info = {}
        self.target = {}

        MF = 0

        # Determine MAT number for this evaluation
        while MF == 0:
            position = fh.tell()
            line = fh.readline()
            MF = int(line[70:72])
        self.material = int(line[66:70])
        fh.seek(position)

        while True:
            # Find next section
            while True:
                position = fh.tell()
                line = fh.readline()
                MAT = int(line[66:70])
                MF = int(line[70:72])
                MT = int(line[72:75])
                if MT > 0 or MAT == 0:
                    fh.seek(position)
                    break

            # If end of material reached, exit loop
            if MAT == 0:
                fh.readline()
                break

            section_data = ''
            while True:
                line = fh.readline()
                if line[72:75] == '  0':
                    break
                else:
                    section_data += line
            self.section[MF, MT] = section_data
        self._read_header()

    def _read_header(self):
        file_obj = StringIO(self.section[1, 451])

        # Information about target/projectile
        items = get_cont_record(file_obj)
        Z, A = divmod(items[0], 1000)
        self.target['atomic_number'] = Z
        self.target['mass_number'] = A
        self.target['mass'] = items[1]
        assert items[2] == 0
        nz = items[3]
        assert items[4] == -1
        ntw = items[5]
        header, items = get_list_record(file_obj)
        self.temperature = header[0]
        assert header[1] == 0.0
        ngn = header[2]
        ngg = header[3]
        nw = header[4]
        assert header[5] == 0
        items = iter(items[ntw:])
        self.sigma_zeros = list(it.islice(items, nz))
        self.neutron_energies = list(it.islice(items, ngn + 1))
        self.photon_energies = list(it.islice(items, ngg + 1))
