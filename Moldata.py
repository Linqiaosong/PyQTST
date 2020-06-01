# encoding=utf-8


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

if __name__ == "__main__":
    A=Moldata(100,200,2)
    B=Moldata(200,100,2)
    C=A+B
    print(C.get_u0k())
    print(C.get_gtk())
    print(C.get_q())
        
        


