# encoding=utf-8

from PyQTST import Moldata,Reaction

# Reacant
R=Moldata(
    U0K=0.0,            # electronic energy+ZPE (internal energy at T=0K) kJ/mol
    GTK=0.0,            # Gibbs free energy at reaction temperature kJ/mol
    Q=1120000000.0      # partition function
    )

# Transition State
TS=Moldata(
    U0K=88.6132504999874,
    GTK=80.92,
    Q=8480000000.0
    )

# Product
P1=Moldata(
    U0K=30.0,
    GTK=0.0,
    Q=1533.0
    )

P2=Moldata(
    U0K=10.0,
    GTK=0.0,
    Q=127865.0
    )

# Reaction A->TS->P
reac=Reaction(
    Nmol=1,             # number of reactant
    molR=R,
    molTS=TS,
    molP=P1+P2,
    Temp=300.0,         # reaction temperature K
    iFreq=-1000.0       # imaginary freq of TS cm-1
    )

# Print Result
reac.printf()

# Print Result to text file (.tst)
reac.print2file(output='1-2.tst')

# Show Reaction Energy image
reac.showimg()
