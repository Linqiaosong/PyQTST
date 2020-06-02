# encoding=utf-8

from PyQTST import run_shermo

# run API, return a Reaction object
result=run_shermo(
    Nmol=2,                     # number of reactant
    Rinput='DA-R.log',          # reactant A freq input for shermo
    R2input='DA-R2.log',        # reactant B freq input for shermo
    TSinput='DA-TS.log',        # TS freq input for shermo
    Pinput='DA-P.log',          # product freq input for shermo
    iFreq=-597.2940,            # imaginary freq of TS (cm-1)
    Temp=373.15                 # reaction temperature (K)
    )


# save output to text file (.tst)
result.print2file(output='api2-1.tst')

# show reaction energy change
result.showimg()
