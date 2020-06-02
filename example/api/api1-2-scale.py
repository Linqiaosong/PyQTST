# encoding=utf-8

from PyQTST import run_shermo

# run API, return a Reaction object
result=run_shermo(
    Nmol=1,                     # number of reactant
    Rinput='EtOH-R.log',        # reactant freq input for shermo
    TSinput='EtOH-TS.log',      # TS freq input for shermo
    Pinput='EtOH-P.log',        # product A freq input for shermo
    P2input='EtOH-P2.log',      # product B freq input for shermo
    iFreq=-2000.2842,           # imaginary freq of TS (cm-1)
    Temp=298.15,                # reaction temperature (K)
    RE=-155.006964,             # reactant high level electronic energy (Eh)
    TSE=-154.8912993,           # TS high level electronic energy (Eh)
    PE=-78.5623346,             # product A high level electronic energy (Eh)
    P2E=-76.4121641,            # product B high level electronic energy (Eh)
    sclZPE=0.9813,              # ZPE freq scale
    sclHeat=1.0004,             # H(T)-H(0) freq scale
    sclS=1.0029                 # S(T) freq scale
    )

# save output to text file (.tst)
result.print2file(
    output='api1-2-scale.tst',  # output filename
    QMethod=False               # not to use partition function method
    )

# show reaction energy change
result.showimg()
