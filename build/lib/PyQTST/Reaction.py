# encoding=utf-8

import sys
import numpy as np
import matplotlib.pyplot as plt
from .Moldata import *


class Reaction:

    nmol=2
    temp=300.0
    ifreq=-0.0
    tunit='k'
    funit='cm-1'

    def __init__(
        self,
        Nmol,
        molR,
        molTS,
        molP,
        Temp=300.0,
        TUnit='K',
        iFreq=-0.1,
        FUnit='cm-1'
        ):

        self.nmol=int(Nmol)
        self.tunit=str.lower(TUnit)
        self.funit=str.lower(FUnit)
        self.molr=molR
        self.molts=molTS
        self.molp=molP

        if self.tunit=='k':
            self.temp=float(Temp)
        elif self.tunit=='c':
            self.temp=c2k(float(Temp))
        elif self.tunit=='f':
            self.temp=f2k(float(Temp))
        else:
            print("The unit of temperature was not defined!")
            exit()

        if self.funit=='cm-1':
            self.ifreq=float(iFreq)
        elif self.funit=='hz':
            self.ifreq=hz2cm(float(iFreq))
        else:
            print("The unit of frequency was not defined!")
            exit()

    def get_nmol(self):
        return self.nmol

    def get_temp(self):
        return self.temp

    def get_ifreq(self):
        return self.ifreq

    def get_ubarrier(self):
        ubarrier=self.molts.get_u0k()-self.molr.get_u0k()
        return ubarrier

    def get_gbarrier(self):
        gbarrier=self.molts.get_gtk()-self.molr.get_gtk()
        return gbarrier

    def get_vfactor(self):
        diffu=self.molp.get_u0k()-self.molr.get_u0k()
        if diffu>=0:
            vfactor=diffu
        else:
            vfactor=0.0
        return vfactor
    
    def get_kappa(self):
        alpha=2.0*pi/h/(-self.ifreq*30000000000.0)
        beta=1.0/(kB*self.temp)
        V=self.get_vfactor()
        barrier=self.get_ubarrier()
        approximate=beta*pi/alpha/np.sin(beta*pi/alpha)
        if alpha>beta:
            kappa=approximate-beta/(alpha-beta)*np.exp((beta-alpha)*(barrier-V)*1000.0/L)
        else:
            kappa=beta/(beta-alpha)*(np.exp((beta-alpha)*(barrier-V)*1000.0/L)-1.0)
        return kappa

    def get_tuneleff(self):
        kappa=self.get_kappa()
        return (kappa-1)/kappa

    def get_kinetic_g(self):
        kappa=self.get_kappa()
        free_energy=self.get_gbarrier()
        if self.nmol==1:
            gk=kappa*(kB*self.temp/h)*np.exp(-free_energy*1000.0/(kB*L*self.temp))
        elif self.nmol==2:
            gk=kappa*kB*self.temp/h*(kB*self.temp/100000.0)*1000000.0*np.exp(-free_energy*1000.0/(kB*L*self.temp))*L*0.001
        return gk

    def get_kinetic_q(self):
        kappa=self.get_kappa()
        barrier=self.get_ubarrier()
        qr=self.molr.get_q()
        qts=self.molts.get_q()
        if self.nmol==1:
            qk=kappa*kB*self.temp/h*qts/qr*np.exp(-barrier*1000.0/(kB*L*self.temp))
        elif self.nmol==2:
            qk=kappa*kB*self.temp/h*(kB*self.temp/100000.0)*1000000.0*qts/qr*np.exp(-barrier*1000.0/(kB*L*self.temp))*L*0.001
        return qk

    def get_halflife_g(self):
        if self.nmol==1:
            gk=self.get_kinetic_g()
            gt=np.log(2)/gk
        elif self.nmol==2:
            gt=0
        return gt

    def get_halflife_q(self):
        if self.nmol==1:
            qk=self.get_kinetic_q()
            qt=np.log(2)/qk
        elif self.nmol==2:
            qt=0
        return qt

    def printf(
        self,
        GMethod=True,
        QMethod=True
        ):

        gr=self.molr.get_gtk()
        gts=self.molts.get_gtk()
        gp=self.molp.get_gtk()
        dg=self.get_gbarrier()
        gk=self.get_kinetic_g()
        gt=self.get_halflife_g()
        qr=self.molr.get_q()
        qts=self.molts.get_q()
        qp=self.molp.get_q()
        qk=self.get_kinetic_q()
        qt=self.get_halflife_q()
        t=self.get_temp()
        n=self.get_nmol()
        ur=self.molr.get_u0k()
        uts=self.molts.get_u0k()
        up=self.molp.get_u0k()
        du=self.get_ubarrier()
        kappa=self.get_kappa()
        eff=self.get_tuneleff()
        dur=self.molp.get_u0k()-self.molr.get_u0k()
        dgr=self.molp.get_gtk()-self.molr.get_gtk()
        kp=np.exp(-dgr/(kB*L*t))

        if GMethod==False and QMethod==True:
            gr=0.0
            gts=0.0
            gp=0.0
            dg=0.0
            gk=0.0
            gt=0.0
        elif GMethod==True and QMethod==False:
            qr=0.0
            qts=0.0
            qp=0.0
            qk=0.0
            qt=0.0
        elif GMethod==False and QMethod==False:
            print('At least one of GMethod and QMethod is True!')
            exit()
        else:
            pass



        print(f'''
=========================================================
                   Calculation Report
                       PyQTST v2.0
            Q. Lin, Wuhan University, 2020
=========================================================


   d8888b. db    db  .d88b.  d888888b .d8888. d888888b 
   88  `8D `8b  d8' .8P  Y8. `~~88~~' 88'  YP `~~88~~' 
   88oodD'  `8bd8'  88    88    88    `8bo.      88    
   88~~~      88    88    88    88      `Y8b.    88    
   88         88    `8P  d8'    88    db   8D    88    
   88         YP     `Y88'Y8    YP    `8888Y'    YP    


PyQTST:

    A software package for calculating chemical reaction
    rate constant by using transition state theory.

Github Website: https://github.com/Linqiaosong/PyQTST

PyPI Website: https://pypi.org/project/PyQTST

Online Documents: https://github.com/Linqiaosong/PyQTST/wiki

PyQTST Citation:

    Q. Lin, PyQTST (version), 
    https://github.com/Linqiaosong/PyQTST, (year)

API for Shermo:

    An API of PyQTST for Shermo software.

Shermo Website: http://sobereva.com/soft/shermo/

Shermo Citation:

    T. Lu, Q. Chen, Shermo: A general code for calculating
    molecular thermodynamic properties, ChemRxiv (2020),
    DOI: 10.26434/chemrxiv.12278801
 



---------------------------------------------------------
                         Part I
             Reaction Molecule Infomations
---------------------------------------------------------

Reactant:

    Electronic Energy + Zero Point Energy:

        U(0 K)=EE+ZPE={ur:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K:

        G({t:.2f} K)={gr:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)/NA={qr:.8E}

Transition State:

    Imaginary Frequency:

        freq={self.ifreq:.2f} cm-1

    Electronic Energy + Zero Point Energy:

        U(0 K)=EE+ZPE={uts:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K without Imaginary
    Frequency:

        G({t:.2f} K)={gts:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)/NA={qts:.8E}

Product:

    Electronic Energy + Zero Point Energy:

        U(0 K)=EE+ZPE={up:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K:

        G({t:.2f} K)={gp:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)/NA={qp:.8E}

---------------------------------------------------------




---------------------------------------------------------
                         Part II
            Reaction Thermodynamics Infomations
---------------------------------------------------------

Pressure Thermodynamic Reference State:

    p*=1.0E+05 Pa

Electronic energy + Single Point Energy Change:

    drUm*(0 K)={dur:.1f} kJ/mol

Gibbs Free Energy Change:

    drGm*({t:.2f} K)={dgr:.1f} kJ/mol

Thermodynamic Equilibrium Constant:

    Kp*({t:.2f} K)={kp:.3E}

---------------------------------------------------------



---------------------------------------------------------
                        Part III
             Reaction Kinetics Infomations
---------------------------------------------------------

Number of Reactive Molecule(s):

    Nmol={n:1d}

Temperature of Reaction:

    T={t:.2f} K

Pressure Thermodynamic Reference State:

    p*=1.0E+05 Pa

Quantum Tunneling Transmission Coefficient:

    Kappa={kappa:.2f}

Contribution of Quantum Tunneling Effect in Reaction:

    eta={100*eff:.2f}%

Reaction Energy Barrier dU(0 K):

    dU={du:.1f} kJ/mol

Reaction Gibbs Free Energy Barrier dG(T={t:.2f} K):

    dG({t:.2f} K)={dg:.1f} kJ/mol        ''')

        if n==1:
            print(f'''
Reaction Rate Constant:

    Gibbs Free Energy Method:

        k({t:.2f} K, GMethod)={gk:.3E} s-1

    Partition Function Method:

        k({t:.2f} K, QMethod)={qk:.3E} s-1

Halflife of Reactant:

    Gibbs Free Energy Method:

        t1/2({t:.2f} K, GMethod)={gt:.3E} s

    Partition Function Method:

        t1/2({t:.2f} K, QMethod)={qt:.3E} s

---------------------------------------------------------
            ''')
        elif n==2:
            print(f'''
Reaction Rate Constant:

    Gibbs Free Energy Method:

        k({t:.2f} K, GMethod)={gk:.3E} (mol/L)-1*s-1

    Partition Function Method:

        k({t:.2f} K, QMethod)={qk:.3E} (mol/L)-1*s-1

---------------------------------------------------------
            ''')
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

    def print2file(
        self,
        output='result.out',
        GMethod=True,
        QMethod=True
        ):

        init = sys.stdout
        sys.stdout=open(output, mode='w')

        self.printf(
            GMethod=GMethod,
            QMethod=QMethod
            )
        
        sys.stdout.close()
        sys.stdout=init

    def showimg(
        self,
        dUimage=True,
        dGimage=True
        ):

        Xaxis=[1,2,2.5,3.5,4,5]
        
        UY1=self.molr.get_u0k()
        GY1=self.molr.get_gtk()
        UY2=self.molts.get_u0k()-UY1
        UY3=self.molp.get_u0k()-UY1
        GY2=self.molts.get_gtk()-GY1
        GY3=self.molp.get_gtk()-GY1
        UY1=0.0
        GY1=0.0
        UYaxis=[UY1,UY1,UY2,UY2,UY3,UY3]
        GYaxis=[GY1,GY1,GY2,GY2,GY3,GY3]

        if dUimage==True and dGimage==False:
            plt.plot(Xaxis,UYaxis)
            plt.title('Electronic Energy + Zero Point Energy')
            plt.ylabel('U(0 K)/(kJ/mol)')
            plt.xticks([])
            plt.show()
        elif dUimage==False and dGimage==True:
            plt.plot(Xaxis,GYaxis)
            plt.title('Gibbs Free Energy')
            plt.ylabel(f'G(T={self.temp:.2f} K)/(kJ/mol)')
            plt.xticks([])
            plt.show()
        elif dUimage==True and dGimage==True:
            plt.figure(1)
            plt.plot(Xaxis,UYaxis)
            plt.title('Electronic Energy + Zero Point Energy')
            plt.ylabel('U(0 K)/(kJ/mol)')
            plt.xticks([])
            plt.figure(2)
            plt.plot(Xaxis,GYaxis)
            plt.title('Gibbs Free Energy')
            plt.ylabel(f'G(T={self.temp:.2f} K)/(kJ/mol)')
            plt.xticks([])
            plt.show()
        result=[Xaxis,UYaxis,GYaxis]
        return result



if __name__ == "__main__":
    
    R=Moldata(
        U0K=0.0,
        GTK=0.0,
        Q=1120000000.0
        )

    R2=Moldata(
        U0K=0.0,
        GTK=0.0,
        Q=80500.0
        )

    TS=Moldata(
        U0K=88.6132504999874,
        GTK=80.92,
        Q=8480000000.0
        )

    P=Moldata(
        U0K=30.0,
        GTK=0.0,
        Q=1.0
        )

    reac=Reaction(
        Nmol=2,
        molR=R+R2,
        molTS=TS,
        molP=P,
        Temp=300.0,
        iFreq=-1000.0
        )

    reac.print2file(GMethod=False)

    reac.showimg()

