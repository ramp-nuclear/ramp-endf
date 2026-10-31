import gc
import warnings

from ramp_endf import Evaluation

# smallest material Evaluation reads: tape ID, MF1/MT451 header (NFOR=6, decay sublibrary), end records
MINIMAL_ENDF = """\
                                                                     1 0  0    0
 1.102400+4 2.379000+1         -1          0          0          01125 1451    1
 0.000000+0 1.000000+0          0          0          0          61125 1451    2
 0.000000+0 0.000000+0          0          0          4          81125 1451    3
 0.000000+0 0.000000+0          0          0          0          01125 1451    4
                                                                  1125 1  099999
                                                                  1125 0  0    0
                                                                     0 0  0    0
                                                                    -1 0  0    0
"""


def test_evaluation_closes_the_file_it_opens(tmp_path):
    path = tmp_path / "minimal.endf"
    path.write_text(MINIMAL_ENDF)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        Evaluation(path)
        gc.collect()
    assert not [w for w in caught if issubclass(w.category, ResourceWarning)]
