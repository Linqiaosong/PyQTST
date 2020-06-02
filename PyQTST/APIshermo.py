# encoding=utf-8

import os
from .Moldata import Moldata
from .Reaction import Reaction

def shermo(input='', Temp=300.0, sclZPE=1.0, sclHeat=1.0, sclS=1.0, E='N/A', shermo_path='Shermo'):
    if E=='N/A':
        command=shermo_path+' '+input+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(sclZPE)+' -sclheat '+str(sclHeat)+' -sclS '+str(sclS)
    else:
        command=shermo_path+' '+input+' -E '+str(E)+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(sclZPE)+' -sclheat '+str(sclHeat)+' -sclS '+str(sclS)
    os.system(command+'> tmp')
    file=open('tmp','r')
    data=file.readlines()[-4:]
    file.close()
    os.remove('tmp')
    data[0]=data[0].split()
    data[3]=data[3].split()
    result={}
    result['U0K']=data[0][11]
    result['GTK']=data[3][9]
    return result
    
def run_shermo(Nmol=1, Rinput='', R2input='N/A', TSinput='', Pinput='', P2input='N/A', iFreq=-0.0, Temp=300.0, RE='N/A', TSE='N/A', PE='N/A', sclZPE=1.0, sclHeat=1.0, sclS=1.0, R2E='N/A', P2E='N/A', shermo_path='Shermo'):
    if Nmol==1:
        Rresult=shermo(input=Rinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=RE, shermo_path=shermo_path)
        TSresult=shermo(input=TSinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=TSE, shermo_path=shermo_path)
        Presult=shermo(input=Pinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=PE, shermo_path=shermo_path)
        R=Moldata(U0K=Rresult['U0K'],GTK=Rresult['GTK'],EUnit='Eh')
        TS=Moldata(U0K=TSresult['U0K'],GTK=TSresult['GTK'],EUnit='Eh')
        if P2input=='N/A':
            P=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
        else:
            P2result=shermo(input=Pinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=PE, shermo_path=shermo_path)
            P1=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
            P2=Moldata(U0K=P2result['U0K'],GTK=P2result['GTK'],EUnit='Eh')
            P=P1+P2            
        cal=Reaction(Nmol=1,molR=R,molTS=TS,molP=P,Temp=Temp,iFreq=iFreq)
        cal.printf(QMethod=False)
    elif Nmol==2 and R2input!='N/A':
        Rresult=shermo(input=Rinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=RE, shermo_path=shermo_path)
        R2result=shermo(input=R2input, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=R2E, shermo_path=shermo_path)
        TSresult=shermo(input=TSinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=TSE, shermo_path=shermo_path)
        Presult=shermo(input=Pinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=PE, shermo_path=shermo_path)
        R=Moldata(U0K=Rresult['U0K'],GTK=Rresult['GTK'],EUnit='Eh')
        R2=Moldata(U0K=R2result['U0K'],GTK=R2result['GTK'],EUnit='Eh')
        TS=Moldata(U0K=TSresult['U0K'],GTK=TSresult['GTK'],EUnit='Eh')
        if P2input=='N/A':
            P=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
        else:
            P2result=shermo(input=Pinput, Temp=Temp, sclZPE=sclZPE, sclHeat=sclHeat, sclS=sclS, E=PE, shermo_path=shermo_path)
            P1=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
            P2=Moldata(U0K=P2result['U0K'],GTK=P2result['GTK'],EUnit='Eh')
            P=P1+P2     
        cal=Reaction(Nmol=2,molR=R,molR2=R2,molTS=TS,molP=P,Temp=Temp,iFreq=iFreq)
        cal.printf(QMethod=False)
    elif Nmol==2 and R2input=='N/A':
        print("Parameter R2input is necessary! ")
        exit()
    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()
    return cal


if __name__ == "__main__":
    print(run_shermo())
    