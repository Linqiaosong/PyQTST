# encoding=utf-8

from PyQTST import run_shermo

def kcal2kj(kcal):
    return float(kcal)*4.1840

def eh2kj(eh):
    return float(eh)*2625.49962

# electronic energy calculated in M052X/6-31G* (Eh) for smd
RgasE=-154.9998078                                            # M052X/6-31G*
RsmdE=-155.0075373                                            # M052X/6-31G* smd
TSgasE=-154.8791114                                           # M052X/6-31G*
TSsmdE=-154.887072                                            # M052X/6-31G* smd
PgasE=-78.5649782                                             # M052X/6-31G*
PsmdE=-78.5631442                                             # M052X/6-31G* smd
P2gasE=-76.3938579                                            # M052X/6-31G*
P2smdE=-76.4081207                                            # M052X/6-31G* smd


# run API, return a Reaction object
result=run_shermo(
    Nmol=1,                                                   # number of reactant
    Rinput='EtOH-R.log',                                      # reactant freq input for shermo
    TSinput='EtOH-TS.log',                                    # TS freq input for shermo
    Pinput='EtOH-P.log',                                      # product A freq input for shermo
    P2input='EtOH-P2.log',                                    # product B freq input for shermo
    iFreq=-2000.2842,                                         # imaginary freq of TS (cm-1)
    Temp=298.15,                                              # reaction temperature (K)
    RE=-155.006964,                                           # reactant high level electronic energy (Eh)
    TSE=-154.8912993,                                         # TS high level electronic energy (Eh)
    PE=-78.5623346,                                           # product A high level electronic energy (Eh)
    P2E=-76.4121641,                                          # product B high level electronic energy (Eh)
    sclZPE=0.9813,                                            # ZPE freq scale
    sclHeat=1.0004,                                           # H(T)-H(0) freq scale
    sclS=1.0029                                               # S(T) freq scale
    )


# add dissolution free energy (kJ/mol)
result.molr.gtk=result.molr.gtk+eh2kj(RsmdE-RgasE)+kcal2kj(1.89)
result.molts.gtk=result.molts.gtk+eh2kj(TSsmdE-TSgasE)+kcal2kj(1.89)
result.molp.gtk=result.molp.gtk+eh2kj(PsmdE-PgasE)+kcal2kj(1.89)+eh2kj(P2smdE-P2gasE)+kcal2kj(1.89)

# print output
result.printf(QMethod=False)    # use Gibbs free energy method only

# save output to text file (.tst)
result.print2file(output='api1-2-smd.tst',QMethod=False)

# show reaction energy change
result.showimg()

