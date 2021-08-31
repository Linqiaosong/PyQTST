# encoding=utf-8

import os
import numpy as np
from Moldata import Moldata
from Reaction import Reaction
import matplotlib.pyplot as plt

def shermo(
    input,
    Temp=300.0,
    sclZPE=1.0,
    sclHeat=1.0,
    sclS=1.0,
    E='N/A',
    shermo_path='Shermo'
    ):

    if E=='N/A':
        command=shermo_path+' '+input+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(sclZPE)+' -sclheat '+str(sclHeat)+' -sclS '+str(sclS)
    else:
        command=shermo_path+' '+input+' -E '+str(E)+' -T '+str(Temp)+' -P 1.0 -sclZPE '+str(sclZPE)+' -sclheat '+str(sclHeat)+' -sclS '+str(sclS)
    
    os.system(command+'> tmp')

    file=open('tmp','r')
    tmp=file.readlines()[-13:]
    file.close()
    os.remove('tmp')

    data=['','','']
    data[0]=tmp[-4].split()
    data[1]=tmp[-1].split()
    data[2]=tmp[0].split()

    result={}
    result['U0K']=data[0][-2]
    result['GTK']=data[1][-2]
    result['Q']=data[2][-1]
    return result
    
def run_shermo(
    Nmol,
    Rinput='',
    R2input='N/A',
    TSinput='',
    Pinput='',
    P2input='N/A',
    iFreq=-0.0,
    Temp=300.0,
    RE='N/A',
    R2E='N/A',
    TSE='N/A',
    PE='N/A',
    P2E='N/A',
    sclZPE=1.0,
    sclHeat=1.0,
    sclS=1.0,
    shermo_path='Shermo'
    ):

    if Nmol==1:
        Rresult=shermo(
            input=Rinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS, E=RE,
            shermo_path=shermo_path
            )

        TSresult=shermo(
            input=TSinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=TSE,
            shermo_path=shermo_path
            )
        
        Presult=shermo(
            input=Pinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=PE,
            shermo_path=shermo_path
            )

        R=Moldata(
            U0K=Rresult['U0K'],
            GTK=Rresult['GTK'],
            Q=Rresult['Q'],
            EUnit='Eh'
            )
        
        TS=Moldata(
            U0K=TSresult['U0K'],
            GTK=TSresult['GTK'],
            Q=TSresult['Q'],
            EUnit='Eh'
            )
        

        if P2input=='N/A':
            P=Moldata(
                U0K=Presult['U0K'],
                GTK=Presult['GTK'],
                Q=Presult['Q'],
                EUnit='Eh'
                )
        else:
            P2result=shermo(
                input=Pinput,
                Temp=Temp,
                sclZPE=sclZPE,
                sclHeat=sclHeat,
                sclS=sclS,
                E=PE,
                shermo_path=shermo_path
                )
            
            P1=Moldata(
                U0K=Presult['U0K'],
                GTK=Presult['GTK'],
                Q=Presult['Q'],
                EUnit='Eh'
                )
            
            P2=Moldata(
                U0K=P2result['U0K'],
                GTK=P2result['GTK'],
                Q=P2result['Q'],
                EUnit='Eh'
                )
            
            P=P1+P2       

             
        cal=Reaction(
            Nmol=1,
            molR=R,
            molTS=TS,
            molP=P,
            Temp=Temp,
            iFreq=iFreq
            )

        cal.printf()
    

    elif Nmol==2 and R2input!='N/A':
        Rresult=shermo(
            input=Rinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=RE,
            shermo_path=shermo_path
            )

        R2result=shermo(
            input=R2input,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=R2E,
            shermo_path=shermo_path
            )
        
        TSresult=shermo(
            input=TSinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=TSE,
            shermo_path=shermo_path
            )
        
        Presult=shermo(
            input=Pinput,
            Temp=Temp,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            E=PE,
            shermo_path=shermo_path
            )

        R=Moldata(
            U0K=Rresult['U0K'],
            GTK=Rresult['GTK'],
            Q=Rresult['Q'],
            EUnit='Eh'
            )

        R2=Moldata(
            U0K=R2result['U0K'],
            GTK=R2result['GTK'],
            Q=R2result['Q'],
            EUnit='Eh'
            )

        TS=Moldata(
            U0K=TSresult['U0K'],
            GTK=TSresult['GTK'],
            Q=TSresult['Q'],
            EUnit='Eh'
            )

        if P2input=='N/A':
            P=Moldata(
                U0K=Presult['U0K'],
                GTK=Presult['GTK'],
                Q=Presult['Q'],
                EUnit='Eh'
                )
        else:
            P2result=shermo(
                input=Pinput,
                Temp=Temp,
                sclZPE=sclZPE,
                sclHeat=sclHeat,
                sclS=sclS,
                E=PE,
                shermo_path=shermo_path
                )
                
            P1=Moldata(
                U0K=Presult['U0K'],
                GTK=Presult['GTK'],
                Q=Presult['Q'],
                EUnit='Eh'
                )

            P2=Moldata(
                U0K=P2result['U0K'],
                GTK=P2result['GTK'],
                Q=P2result['Q'],
                EUnit='Eh'
                )
            
            P=P1+P2     

        
        cal=Reaction(
            Nmol=2,
            molR=R+R2,
            molTS=TS,
            molP=P,
            Temp=Temp,
            iFreq=iFreq
            )

        cal.printf()

    
    elif Nmol==2 and R2input=='N/A':
        print("Parameter R2input is necessary! ")
        exit()
    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()
    return cal


def scan_shermo(
    Nmol,
    Rinput='',
    R2input='N/A',
    TSinput='',
    Pinput='',
    P2input='N/A',
    iFreq=-0.0,
    Temp=[],
    RE='N/A',
    R2E='N/A',
    TSE='N/A',
    PE='N/A',
    P2E='N/A',
    sclZPE=1.0,
    sclHeat=1.0,
    sclS=1.0,
    shermo_path='Shermo'
    ):

    result=[]
    Xaxis=Temp
    GYaxis=[]
    QYaxis=[]

    for temp in Temp:
        calc=run_shermo(
            Nmol=Nmol,
            Rinput=Rinput,
            R2input=R2input,
            TSinput=TSinput,
            Pinput=Pinput,
            P2input=P2input,
            iFreq=iFreq,
            Temp=float(temp),
            RE=RE,
            R2E=R2E,
            TSE=TSE,
            PE=PE,
            P2E=P2E,
            sclZPE=sclZPE,
            sclHeat=sclHeat,
            sclS=sclS,
            shermo_path=shermo_path
            )
        GYaxis.append(calc.get_kinetic_g())
        QYaxis.append(calc.get_kinetic_q())
    
    xx=np.linspace(min(Temp),max(Temp),200)
    GFunc=np.poly1d(np.polyfit(Xaxis,GYaxis,3))
    QFunc=np.poly1d(np.polyfit(Xaxis,QYaxis,3))

    result.append(Temp)
    result.append(GYaxis)
    result.append(QYaxis)

    if Nmol==1:
        print(f'''

---------------------------------------------------------
                        Part IV
  Relationship between Reaction Kinetics and Temperature
---------------------------------------------------------

Gibbs Free Energy Method:

    k(T)/s-1
    
    =  ({GFunc.coeffs[0]:.3E})\t * (T/K)^3
    
     + ({GFunc.coeffs[1]:.3E})\t * (T/K)^2
    
     + ({GFunc.coeffs[2]:.3E})\t * (T/K)
    
     + ({GFunc.coeffs[3]:.3E})


    RMSD=


Partition Function Method:

    k(T)/s-1
    
    =  ({QFunc.coeffs[0]:.3E})\t * (T/K)^3
    
     + ({QFunc.coeffs[1]:.3E})\t * (T/K)^2
    
     + ({QFunc.coeffs[2]:.3E})\t * (T/K)
    
     + ({QFunc.coeffs[3]:.3E})


    RMSD=

---------------------------------------------------------
        ''')

    elif Nmol==2:
        print(f'''

---------------------------------------------------------
                        Part IV
  Relationship between Reaction Kinetics and Temperature
---------------------------------------------------------

Gibbs Free Energy Method:

    k(T)/((mol/L)-1*s-1)

    =  ({GFunc.coeffs[0]:.3E})\t * (T/K)^3

     + ({GFunc.coeffs[1]:.3E})\t * (T/K)^2

     + ({GFunc.coeffs[2]:.3E})\t * (T/K)

     + ({GFunc.coeffs[3]:.3E})


    RMSD=


Partition Function Method:

    k(T)/((mol/L)-1*s-1)

    =  ({QFunc.coeffs[0]:.3E})\t * (T/K)^3

     + ({QFunc.coeffs[1]:.3E})\t * (T/K)^2

     + ({QFunc.coeffs[2]:.3E})\t * (T/K)

     + ({QFunc.coeffs[3]:.3E})


    RMSD=

---------------------------------------------------------
        ''')

    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()


    print('''

d8b   db  .d88b.  d8888b. .88b  d88.  .d8b.  db      
888o  88 .8P  Y8. 88  `8D 88'YbdP`88 d8' `8b 88      
88V8o 88 88    88 88oobY' 88  88  88 88ooo88 88      
88 V8o88 88    88 88`8b   88  88  88 88~~~88 88      
88  V888 `8b  d8' 88 `88. 88  88  88 88   88 88booo. 
VP   V8P  `Y88P'  88   YD YP  YP  YP YP   YP Y88888P 
                                                     
                                                     
d88888b d8b   db d8888b.      db                     
88'     888o  88 88  `8D      88                     
88ooooo 88V8o 88 88   88      YP                     
88~~~~~ 88 V8o88 88   88                             
88.     88  V888 88  .8D      db                     
Y88888P VP   V8P Y8888D'      YP                   

            ''')

    plt.plot(Xaxis,GYaxis,'rx',xx,GFunc(xx),'-')
    plt.title('Reaction Rate Constant (Gibbs Free Energy Method)')
    if Nmol==1:
        plt.ylabel('k/(s-1)')
    elif Nmol==2:
        plt.ylabel('k/((mol/L)-1*s-1)')
    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()
    plt.xlabel('Temperature/K')
    plt.figure(2)
    plt.plot(Xaxis,QYaxis,'rx',xx,QFunc(xx),'-')
    plt.title('Reaction Rate Constant (Partition Function Method)')
    if Nmol==1:
        plt.ylabel('k/(s-1)')
    elif Nmol==2:
        plt.ylabel('k/((mol/L)-1*s-1)')
    else:
        print("The number of molecules in primitive chemical reaction is not suitable! ")
        exit()
    plt.xlabel('Temperature/K')
    plt.show()


    return result




if __name__ == "__main__":
    result=scan_shermo(
        Nmol=2,                     
        Rinput='DA-R.log',   
        R2input='DA-R2.log',
        TSinput='DA-TS.log',      
        Pinput='DA-P.log',          
        iFreq=-2000.2842,           
        Temp=list(range(273,373))
        )
    print(result)
    
