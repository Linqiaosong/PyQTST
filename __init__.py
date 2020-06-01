# encoding=utf-8

from Moldata import Moldata
from Reaction import Reaction
from APIshermo import APIshermo

R=Moldata(U0K=0.0,GTK=0.0,Q=1120000000.0)
R2=Moldata(U0K=0.0,GTK=0.0,Q=80500.0)
TS=Moldata(U0K=88.6132504999874,GTK=80.92,Q=8480000000.0)
P=Moldata(U0K=30.0,GTK=0.0,Q=1.0)
reac=Reaction(Nmol=2,molR=R,molR2=R2,molTS=TS,molP=P,Temp=300.0,iFreq=-1000.0)
reac.printf()
reac.showimg()