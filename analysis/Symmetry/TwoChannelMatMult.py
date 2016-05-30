from sympy import *
from sympy.matrices import *

def getJostSymbols(m,n):
    return symbols('k'+str(m)+' '+'k'+str(n)+' '+'A'+str(m)+str(n)+' '+'B'+str(m)+str(n))

def inJost(m,n):
    km,kn,Amn,Bmn = getJostSymbols(m,n)
    return kn/km*Amn - kn*Bmn*1.0j

def outJost(m,n):
    km,kn,Amn,Bmn = getJostSymbols(m,n)
    return kn/km*Amn + kn*Bmn*1.0j

def getNum():
    return Matrix([[outJost(1,1),outJost(1,2)], [outJost(2,1),outJost(2,2)]])

def getDenum():
    return Matrix([[inJost(1,1),inJost(1,2)], [inJost(2,1),inJost(2,2)]])


def getPPDict():
    val = 1.0
    k1val = 1.0
    k2val = 3.0
    return {'k1': k1val, 'k2': k2val, 'A11': val, 'A12': val, 'A21': val, 'A22': val, 'B11': val, 'B12': val, 'B21': val, 'B22': val }

def getMMDict():
    val = 1.0
    k1val = 1.0
    k2val = 3.0
    return {'k1': -k1val, 'k2': -k2val, 'A11': val, 'A12': val, 'A21': val, 'A22': val, 'B11': val, 'B12': val, 'B21': val, 'B22': val }

#S = getNum() * getDenum().inv()
#print(latex(S[1,0]))

S = getNum() * getDenum().inv()
print S[0,0].evalf(subs=getPPDict())
print S[0,0].evalf(subs=getMMDict())
print "\n"
print S[0,1].evalf(subs=getPPDict())
print S[0,1].evalf(subs=getMMDict())
print "\n"
print S[1,0].evalf(subs=getPPDict())
print S[1,0].evalf(subs=getMMDict())
print "\n"
print S[1,1].evalf(subs=getPPDict())
print S[1,1].evalf(subs=getMMDict())

#S = Matrix( [[Spp[0,0].as_real_imag(),Spp[0,1].as_real_imag()], [Spp[1,0].as_real_imag(),Spp[1,1].as_real_imag()]] )

#print(latex(S))
#print(pretty(res))