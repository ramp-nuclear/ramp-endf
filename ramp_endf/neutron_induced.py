"""Data used to understand ENDF cross sections (neutron based)

"""
from enum import IntEnum


class InducedReactionType(IntEnum):
    """ENDF reaction numbers.
    see mcnp manual G-2, or ENDF6 formats manual appendix B

    We currently support up to MT=250
    """
    TOTAL = 1
    ELASTIC = 2
    NON_ELASTIC = 3
    NN = 4  # Production of one neutron in the exit channel. Sum of 50-91
    OTHER = 5  # Sum of all others not in an MT. Needed for the sum to equal MT=1.
    CONTINUUM = 10  # Total continuum reactions, excluding discretes.
    ND2N = 11  # Production of two neutrons and a deutron, plus residual
    N2N = 16  # (n,2n) cross section
    N3N = 17  # (n,3n) cross section
    FISS = 18  # fission cross section
    NF = 19  # fission cross section? (n,f) in the ENDF manual page B.3
    SECOND_FISSION = 20  # Second-chance fission (not defined for charged particles)
    THIRD_FISSION = 21  # Third-chance fission (not defined for charged particles)
    NNALPHA = 22  # Neutron and an Alpha, plus residual
    NN3ALPHA = 23  # Production of a neutron and 3 alpha particles, plus residual
    N2NALPHA = 24  # Production of 2 neutrons and an alpha particle, plus residual
    N3NALPHA = 25  # Production of 3 neutrons and an alpha particle, plus residual
    OBSOLETE_ISOMERIC = 26  # Used in ENDF version 5 only for isomeric state in files 6, 8, 9, 10.
    NABS = 27  # Sum of 18 and 102-117, rarely used.
    NNP = 28  # Production of a proton and a neutron, plus residual
    NN2ALPHA = 29  # Production of a neutron and 2 alpha particles, plus residual
    N2N2ALPHA = 30  # Production of 2 neutrons and 3 alpha particles, plus residual
    NDN = 32  # Production of a neutron and a deutron, plus residual
    NTN = 33  # Production of a neutron and a triton, plus residual
    NNHe3 = 34  # Production of a neutron and a He3, plus residual
    NND2ALPHA = 35  # Production of a neutron, a deutron and 2 alpha particles, plus residual
    NNT2ALPHA = 36  # Production of a neutron, a triton and 2 alpha particles, plus residual
    N4N = 37  # Production of 4 neutrons
    FOURTH_FISSION = 38  # Fourth-chance fission (not defined for charged particles)
    N2NP = 41  # Two neutrons and a proton
    N3NP = 42  # Three neutrons and a proton
    NN2P = 44  # A neutron and 2 protons
    NNPALPHA = 45  # A neutron, a proton and an alpha particle
    YN_GROUND = 50  # A neutron, leaving the residual in the ground state. Not allowed for incident neutrons, use MT=2
    NN_1 = 51  # A neutron, leaving the residual in the 1 excited state.
    NN_2 = 52  # A neutron, leaving the residual in the 2 excited state.
    NN_3 = 53  # A neutron, leaving the residual in the 3 excited state.
    NN_4 = 54  # A neutron, leaving the residual in the 4 excited state.
    NN_5 = 55  # A neutron, leaving the residual in the 5 excited state.
    NN_6 = 56  # A neutron, leaving the residual in the 6 excited state.
    NN_7 = 57  # A neutron, leaving the residual in the 7 excited state.
    NN_8 = 58  # A neutron, leaving the residual in the 8 excited state.
    NN_9 = 59  # A neutron, leaving the residual in the 9 excited state.
    NN_10 = 60  # A neutron, leaving the residual in the 10 excited state.
    NN_11 = 61  # A neutron, leaving the residual in the 11 excited state.
    NN_12 = 62  # A neutron, leaving the residual in the 12 excited state.
    NN_13 = 63  # A neutron, leaving the residual in the 13 excited state.
    NN_14 = 64  # A neutron, leaving the residual in the 14 excited state.
    NN_15 = 65  # A neutron, leaving the residual in the 15 excited state.
    NN_16 = 66  # A neutron, leaving the residual in the 16 excited state.
    NN_17 = 67  # A neutron, leaving the residual in the 17 excited state.
    NN_18 = 68  # A neutron, leaving the residual in the 18 excited state.
    NN_19 = 69  # A neutron, leaving the residual in the 19 excited state.
    NN_20 = 70  # A neutron, leaving the residual in the 20 excited state.
    NN_21 = 71  # A neutron, leaving the residual in the 21 excited state.
    NN_22 = 72  # A neutron, leaving the residual in the 22 excited state.
    NN_23 = 73  # A neutron, leaving the residual in the 23 excited state.
    NN_24 = 74  # A neutron, leaving the residual in the 24 excited state.
    NN_25 = 75  # A neutron, leaving the residual in the 25 excited state.
    NN_26 = 76  # A neutron, leaving the residual in the 26 excited state.
    NN_27 = 77  # A neutron, leaving the residual in the 27 excited state.
    NN_28 = 78  # A neutron, leaving the residual in the 28 excited state.
    NN_29 = 79  # A neutron, leaving the residual in the 29 excited state.
    NN_30 = 80  # A neutron, leaving the residual in the 30 excited state.
    NN_31 = 81  # A neutron, leaving the residual in the 31 excited state.
    NN_32 = 82  # A neutron, leaving the residual in the 32 excited state.
    NN_33 = 83  # A neutron, leaving the residual in the 33 excited state.
    NN_34 = 84  # A neutron, leaving the residual in the 34 excited state.
    NN_35 = 85  # A neutron, leaving the residual in the 35 excited state.
    NN_36 = 86  # A neutron, leaving the residual in the 36 excited state.
    NN_37 = 87  # A neutron, leaving the residual in the 37 excited state.
    NN_38 = 88  # A neutron, leaving the residual in the 38 excited state.
    NN_39 = 89  # A neutron, leaving the residual in the 39 excited state.
    NN_40 = 90  # A neutron, leaving the residual in the 40 excited state.
    NN_CONTINUUM = 91  # A neutron, leaving the residual in the continuum state.
    NDISAPPEAR = 101  # Neutron disappearance. Sum of 102-117
    NG = 102  # (n,gamma) cross section
    NP = 103  # (n,p) cross section
    ND = 104  # (n,d) cross section
    NT = 105  # (n,t) cross section
    NHE3 = 106  # (n,He3) cross section
    NALPHA = 107  # (n,alpha) cross section
    N2Alpha = 108  # (n,2alpha)
    N3Alpha = 109  # (n,3alpha)
    N2P = 111  # (n,2p)
    NPALPHA = 112  # Production of a proton and an alpha particle
    NT2ALPHA = 113  # Production of a triton and 2 alpha particles
    ND2ALPHA = 114  # Production of a deutron and 2 alpha particles
    NPD = 115  # Production of a proton and a deutron
    NPT = 116  # Production of a proton and a triton
    NDALPHA = 117  # Production of a deutron and an alpha particle
    OBSOLETE_DESTRUCTION = 120  # Version 5 only: Nonelastic minus total
    NRES = 151  # Resonance parameters used to calculate cross sections at different temperatures in neutrons only!

    Ntot = 201  # Total neutron production, redundant for derived files only
    NGtot = 202  # Total Gamma production, redundant for derived files only
    NPtot = 203  # Total proton production, redundant for derived files only
    NDtot = 204  # Total deuterium production, redundant for derived files only
    NTtot = 205  # Total tritium production, redundant for derived files only
    NHE3tot = 206  # Total He3 production, redundant for derived files only
    NALPHAtot = 207  # Total Alpha production, redundant for derived files only
    NPiPlustot = 208  # Total Pi+ production, for use in high energy
    NPi0tot = 209  # Total Pi0 production, for use in high energy
    NPiMinustot = 210  # Total Pi- production, for use in high energy
    NMuPlustot = 211  # Total Mu+ production, for use in high energy
    NMuMinustot = 212  # Total Mu- production, for use in high energy
    NKPlustot = 213  # Total K+ production, for use in high energy
    NK0Longtot = 214  # Total K0 (long) production, for use in high energy
    NK0Shorttot = 215  # Total K0 (short) production, for use in high energy
    NKMinustot = 216  # Total K- production, for use in high energy
    NAntiP = 217  # Total anti proton production, for use in high energy
    NAntin = 218  # Total anti neutron production, for use in high energy


