# encoding=utf-8

from .Moldata import Moldata
from .Reaction import Reaction
from .APIshermo import shermo,run_shermo

<<<<<<< HEAD
=======
c=299792458
kB=1.3806503E-23
h=6.6260696E-34
L=6.02214179E+23
R=kB*L
pi=3.14159265358979

def kcal2kj(kcal):
    return float(kcal)*4.1840

def ev2kj(ev):
    return float(ev)*96.485

def eh2kj(eh):
    return float(eh)*2625.49962

def c2k(oc):
    return float(oc)+273.15

def f2k(of):
    return (float(of)-32)/1.8+273.15

def hz2cm(hz):
    return hz/(100*c)

class Moldata:

    u0k=0.0
    gtk=0.0
    q=1.0
    eunit='kj/mol'

    def __init__(self, U0K=0.0, GTK=0.0, Q=1.0, EUnit='kj/mol'):
        self.eunit=str.lower(EUnit)
        self.q=float(Q)
        if self.eunit=='kj/mol':
            self.u0k=float(U0K)
            self.gtk=float(GTK)
        elif self.eunit=='kcal/mol':
            self.u0k=kcal2kj(float(U0K))
            self.gtk=kcal2kj(float(GTK))
        elif self.eunit=='ev':
            self.u0k=ev2kj(float(U0K))
            self.gtk=ev2kj(float(GTK))
        elif self.eunit=='eh':
            self.u0k=eh2kj(float(U0K))
            self.gtk=eh2kj(float(GTK))
        else:
            print("The unit of energy was not defined!")
            exit()

    def __add__(self,other):
        result=Moldata(U0K=self.u0k+other.u0k,GTK=self.gtk+other.gtk,Q=self.q*other.q,EUnit='kj/mol')
        return result

    def get_u0k(self):
        return self.u0k

    def get_gtk(self):
        return self.gtk

    def get_q(self):
        return self.q


class Reaction:

    nmol=2
    temp=300.0
    ifreq=-0.0
    tunit='k'
    funit='cm-1'

    def __init__(self, Nmol, molR, molTS, molP, molR2='N/A', Temp=300.0, TUnit='K', iFreq=-0.0, FUnit='cm-1'):
        self.nmol=int(Nmol)
        self.tunit=str.lower(TUnit)
        self.funit=str.lower(FUnit)
        if self.nmol==1:
            self.molr=molR
        elif self.nmol==2 and molR2!='N/A':
            self.molr=molR
            self.molr2=molR2
        else:
            print("The number of molecules in primitive chemical reaction is not suitable! ")
            exit()
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
        if self.nmol==1:
            ubarrier=self.molts.get_u0k()-self.molr.get_u0k()
        elif self.nmol==2:
            ubarrier=self.molts.get_u0k()-self.molr.get_u0k()-self.molr2.get_u0k()
        return ubarrier

    def get_gbarrier(self):
        if self.nmol==1:
            gbarrier=self.molts.get_gtk()-self.molr.get_gtk()
        elif self.nmol==2:
            gbarrier=self.molts.get_gtk()-self.molr.get_gtk()-self.molr2.get_gtk()
        return gbarrier

    def get_vfactor(self):
        if self.nmol==1:
            diffu=self.molp.get_u0k()-self.molr.get_u0k()
        elif self.nmol==2:
            diffu=self.molp.get_u0k()-self.molr.get_u0k()-self.molr2.get_u0k()
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
        barrier=barrier=self.get_ubarrier()
        qr=self.molr.get_q()
        qts=self.molts.get_q()
        if self.nmol==1:
            qk=kappa*kB*self.temp/h*qts/qr*np.exp(-barrier*1000.0/(kB*L*self.temp))
        elif self.nmol==2:
            qr2=self.molr2.get_q()
            qk=kappa*kB*self.temp/h*(kB*self.temp/100000.0)*1000000.0*qts/qr/qr2*np.exp(-barrier*1000.0/(kB*L*self.temp))*L*0.001
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

    def printf(self,GMethod=True,QMethod=True):

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
        if n==1:
            dgr=self.molp.get_gtk()-self.molr.get_gtk()
        elif n==2:
            dgr=self.molp.get_gtk()-self.molr.get_gtk()-self.molr2.get_gtk()
            gr2=self.molr2.get_gtk()
            qr2=self.molr2.get_q()
            ur2=self.molr2.get_u0k()
        kp=np.exp(-dgr/(kB*L*t))

        if GMethod==False and QMethod==True:
            gr=0.0
            gr2=0.0
            gts=0.0
            gp=0.0
            dg=0.0
            gk=0.0
            gt=0.0
        elif QMethod==False and GMethod==True:
            qr=0.0
            qr2=0.0
            qts=0.0
            qp=0.0
            qk=0.0
            qt=0.0
        elif GMethod==False and QMethod==False:
            print('At least one of GMethod and QMethod is True!')
            exit()
        else:
            pass



        print('''

=========================================================
                   Calculation Report
                       PyQTST v2.0
            Q. Lin, Wuhan University, 2020
=========================================================


        ''')


        if n==1:
            print(f'''


---------------------------------------------------------
                         Part I
             Reaction Molecule Infomations
---------------------------------------------------------

1. Reactant A:

(1) Electronic Energy and Zero Point Energy:

    U(0 K)=EE+ZPE={ur:.1f} kJ/mol

(2) Gibbs Free Energy at T={t:.2f} K:

    G(T={t:.2f} K)={gr:.1f} kJ/mol

(3) Total Partition Function without Zero Point Energy:

    Q(V=0)={qr:.8E}''')
        elif n==2:
            print(f'''


---------------------------------------------------------
                         Part I
             Reaction Molecule Infomations
---------------------------------------------------------

Reactant A:

    Electronic Energy and Zero Point Energy:

        U(0 K)=EE+ZPE={ur:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K:

        G(T={t:.2f} K)={gr:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)={qr:.8E}

Reactant B:

    Electronic Energy and Zero Point Energy:

        U(0 K)=EE+ZPE={ur2:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K:

        G(T={t:.2f} K)={gr2:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)={qr2:.8E}''')
        
        print(f'''
Transition State:

    Electronic Energy and Zero Point Energy:

        U(0 K)=EE+ZPE={uts:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K without Imaginary Frequency:

        G(T={t:.2f} K)={gts:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)={qts:.8E}

Product:

    Electronic Energy and Zero Point Energy:

        U(0 K)=EE+ZPE={up:.1f} kJ/mol

    Gibbs Free Energy at T={t:.2f} K:

        G(T={t:.2f} K)={gp:.1f} kJ/mol

    Total Partition Function without Zero Point Energy:

        Q(V=0)={qp:.8E}

---------------------------------------------------------




---------------------------------------------------------
                         Part II
            Reaction Thermodynamics Infomations
---------------------------------------------------------

Pressure Thermodynamic Reference State:

    p*=1.0E+5 Pa

Gibbs Free Energy Change:

    drGm*(T={t:.2f} K)={dgr:.1f} kJ/mol

Thermodynamic Equilibrium Constant:

    Kp*(T={t:.2f} K)={kp:.3E}

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

    dG(T={t:.2f} K)={dg:.1f} kJ/mol        ''')

        if n==1:
            print(f'''
Reaction Rate Constant:

    Gibbs Free Energy Method:

        k(T={t:.2f} K, GMethod)={gk:.3E} s-1

    Partition Function Method:

        k(T={t:.2f} K, QMethod)={qk:.3E} s-1

Halflife of Reactant:

    Gibbs Free Energy Method:

        t1/2(T={t:.2f} K, GMethod)={gt:.3E} s

    Partition Function Method:

        t1/2(T={t:.2f} K, QMethod)={qt:.3E} s

---------------------------------------------------------
            ''')
        elif n==2:
            print(f'''
Reaction Rate Constant:

    Gibbs Free Energy Method:

        k(T={t:.2f} K, GMethod)={gk:.3E} (mol/L)-1*s-1

    Partition Function Method:

        k(T={t:.2f} K, QMethod)={qk:.3E} (mol/L)-1*s-1

---------------------------------------------------------




Normal End!
            ''')

    def showimg(self,dUimage=True,dGimage=True):
        Xaxis=[1,2,2.5,3.5,4,5]
        
        if self.nmol==1:
            UY1=self.molr.get_u0k()
            GY1=self.molr.get_gtk()
        elif self.nmol==2:
            UY1=self.molr.get_u0k()+self.molr2.get_u0k()
            GY1=self.molr.get_gtk()+self.molr2.get_u0k()
        UY2=self.molts.get_u0k()-UY1
        UY3=self.molp.get_u0k()-UY1
        GY2=self.molts.get_gtk()-GY1
        GY3=self.molp.get_gtk()-GY1
        UY1=0.0
        GY1=0.0
        UYaxis=[UY1,UY1,UY2,UY2,UY3,UY3]
        GYaxis=[GY1,GY1,GY2,GY2,GY3,GY3]
        if dUimage==False and dGimage==True:
            plt.plot(Xaxis,GYaxis)
            plt.title('Electronic Energy + Zero Point Energy')
            plt.ylabel('U(0 K)/(kJ/mol)')
            plt.xticks([])
            plt.show()
        elif dUimage==True and dGimage==False:
            plt.plot(Xaxis,UYaxis)
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
>>>>>>> ab336a70a7dfcd6a517e22f5db707b0a040cf726



if __name__ == "__main__":
    R=Moldata(U0K=0.0,GTK=0.0,Q=1120000000.0)
    R2=Moldata(U0K=0.0,GTK=0.0,Q=80500.0)
    TS=Moldata(U0K=88.6132504999874,GTK=80.92,Q=8480000000.0)
    P=Moldata(U0K=30.0,GTK=0.0,Q=1.0)
    reac=Reaction(Nmol=2,molR=R,molR2=R2,molTS=TS,molP=P,Temp=300.0,iFreq=-1000.0)
    reac.printf()
    reac.showimg()
<<<<<<< HEAD

    
=======
>>>>>>> ab336a70a7dfcd6a517e22f5db707b0a040cf726
