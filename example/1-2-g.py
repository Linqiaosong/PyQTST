# encoding=utf-8

from PyQTST import Moldata,Reaction

# Reacant
R=Moldata(
    GTK=0.0            # Gibbs free energy at reaction temperature kJ/mol
    )

# Transition State
TS=Moldata(
    GTK=80.92
    )

# Product
P1=Moldata(
    GTK=0.0
    )

P2=Moldata(
    GTK=0.0
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
reac.printf(QMethod=False)  # use Gibbs free energy method moly

# Print Result to text file (.tst)
reac.print2file(output='1-2-g.tst', QMethod=False)     # use Gibbs free energy method moly

# Show Reaction Energy image
reac.showimg(dUimage=False)