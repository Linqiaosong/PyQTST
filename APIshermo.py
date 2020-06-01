# encoding=utf-8

import os
from Reaction import Reaction
from Moldata import Moldata

def shermo(output='C2H5_optfreq.out', Temp=300.0, ZPEscale=1.0, Hscale=1.0, Sscale=1.0, shermo_path='Shermo', EE='N/A'):
    if EE=='N/A':
        command=shermo_path+' '+output+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(ZPEscale)+' -sclheat '+str(Hscale)+' -sclS '+str(Sscale)
    else:
        command=shermo_path+' '+output+' -E '+str(EE)+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(ZPEscale)+' -sclheat '+str(Hscale)+' -sclS '+str(Sscale)
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
    
def run_shermo(Nmol=1, Routput='C2H5_optfreq.out', TSoutput='C2H5_optfreq.out', Poutput='C2H5_optfreq.out', imgFreq=-0.0, T=300.0, shermo_Path='Shermo', RE='N/A', TSE='N/A', PE='N/A', sclZPE=1.0, sclHeat=1.0, sclS=1.0, R2output='N/A', R2E='N/A', P2output='N/A', P2E='N/A'):
    if Nmol==1:
        Rresult=shermo(output=Routput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=RE, shermo_path=shermo_Path)
        TSresult=shermo(output=TSoutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=TSE, shermo_path=shermo_Path)
        Presult=shermo(output=Poutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=PE, shermo_path=shermo_Path)
        R=Moldata(U0K=Rresult['U0K'],GTK=Rresult['GTK'],EUnit='Eh')
        TS=Moldata(U0K=TSresult['U0K'],GTK=TSresult['GTK'],EUnit='Eh')
        if P2output=='N/A':
            P=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
        else:
            P2result=shermo(output=Poutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=PE, shermo_path=shermo_Path)
            P1=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
            P2=Moldata(U0K=P2result['U0K'],GTK=P2result['GTK'],EUnit='Eh')
            P=P1+P2            
        cal=Reaction(Nmol=1,molR=R,molTS=TS,molP=P,Temp=T,iFreq=imgFreq)
        cal.printf(QMethod=False)
        cal.showimg()
    elif Nmol==2 and R2output!='N/A':
        Rresult=shermo(output=Routput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=RE, shermo_path=shermo_Path)
        R2result=shermo(output=R2output, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=R2E, shermo_path=shermo_Path)
        TSresult=shermo(output=TSoutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=TSE, shermo_path=shermo_Path)
        Presult=shermo(output=Poutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=PE, shermo_path=shermo_Path)
        R=Moldata(U0K=Rresult['U0K'],GTK=Rresult['GTK'],EUnit='Eh')
        R2=Moldata(U0K=R2result['U0K'],GTK=R2result['GTK'],EUnit='Eh')
        TS=Moldata(U0K=TSresult['U0K'],GTK=TSresult['GTK'],EUnit='Eh')
        if P2output=='N/A':
            P=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
        else:
            P2result=shermo(output=Poutput, Temp=T, ZPEscale=sclZPE, Hscale=sclHeat, Sscale=sclS, EE=PE, shermo_path=shermo_Path)
            P1=Moldata(U0K=Presult['U0K'],GTK=Presult['GTK'],EUnit='Eh')
            P2=Moldata(U0K=P2result['U0K'],GTK=P2result['GTK'],EUnit='Eh')
            P=P1+P2     
        cal=Reaction(Nmol=2,molR=R,molR2=R2,molTS=TS,molP=P,Temp=T,iFreq=imgFreq)
        cal.printf(QMethod=False)
        cal.showimg()
    elif Nmol==2 and R2output=='N/A':
        print("Parameter R2output is necessary! ")
        exit()
    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()

if __name__ == "__main__":
    run_shermo(Nmol=1, Routput='DA-R.log', TSoutput='DA-TS.log', imgFreq=-597.2940, T=300.0, shermo_Path=r'.\Shermo.exe', Poutput='DA-P.log', RE=-234.31053627, TSE=-234.27767276, PE=-234.37295465, sclZPE=1.0, sclHeat=1.0, sclS=1.0)
