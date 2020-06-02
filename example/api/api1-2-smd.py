# encoding=utf-8

from PyQTST import run_shermo

# electronic energy calculated in M052X/6-31G* (Eh)
RgasE=-154.9998078                                            # gas
RsmdE=-155.0075373                                            # smd model
TSgasE=-154.8791114                                           # gas
TSsmdE=-154.887072                                            # smd model
PgasE=-78.5649782                                             # gas
PsmdE=-78.5631442                                             # smd model
P2gasE=-76.3938579                                            # gas
P2smdE=-76.4081207                                            # smd model


# run API, return a Reaction object
result=run_shermo(
    Nmol=1,                                                   # number of reactant
    Rinput='EtOH-R.log',                                      # reactant freq input for shermo
    TSinput='EtOH-TS.log',                                    # TS freq input for shermo
    Pinput='EtOH-P.log',                                      # product A freq input for shermo
    P2input='EtOH-P2.log',                                    # product B freq input for shermo
    iFreq=-2000.2842,                                         # imaginary freq of TS (cm-1)
    Temp=298.15,                                              # reaction temperature (K)
    RE=-155.006964+(RsmdE-RgasE)+(1.89/627.51),               # reactant high level electronic energy (Eh)
    TSE=-154.8912993+(TSsmdE-TSgasE)+(1.89/627.51),           # TS high level electronic energy (Eh)
    PE=-78.5623346+(PsmdE-PgasE)+(1.89/627.51),               # product A high level electronic energy (Eh)
    P2E=-76.4121641+(P2smdE-P2gasE)+(1.89/627.51),            # product B high level electronic energy (Eh)
    sclZPE=0.9813,                                            # ZPE freq scale
    sclHeat=1.0004,                                           # H(T)-H(0) freq scale
    sclS=1.0029                                               # S(T) freq scale
    )

# save output to text file (.tst)
result.print2file(output='api1-2-smd.tst')



# show reaction energy change
result.showimg()
