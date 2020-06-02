# encoding=utf-8

from PyQTST import Moldata,Reaction

# Reacant
R=Moldata(
    U0K=0.0,            # electronic energy+ZPE (internal energy at T=0K) kJ/mol
    GTK=0.0             # Gibbs free energy at reaction temperature kJ/mol
    )

R2=Moldata(
    U0K=0.0,            # electronic energy+ZPE (internal energy at T=0K) kJ/mol
    GTK=0.0             # Gibbs free energy at reaction temperature kJ/mol
    )

# Transition State
TS=Moldata(
    U0K=88.6132504999874,
    GTK=112.0
    )

# Product
P=Moldata(
    U0K=30.0,
    GTK=10.0
    )


# Reaction A->TS->P
reac=Reaction(
    Nmol=2,             # number of reactant
    molR=R+R2,
    molTS=TS,
    molP=P,
    Temp=300.0,         # reaction temperature K
    iFreq=-1000.0       # imaginary freq of TS cm-1
    )


# Print Result
reac.printf(QMethod=False) # use Gibbs free energy method moly

# Print Result to text file (.tst)
reac.print2file(output='2-1-g.tst',QMethod=False) # use Gibbs free energy method moly

# Show Reaction Energy image
reac.showimg()
