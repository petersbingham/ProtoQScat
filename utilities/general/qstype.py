import cmath
import sympy.mpmath as mpmath
import numpy as np
import sympy as sy

MODE_NORM = 0
MODE_MPMATH = 1

##########################################################
################### Configuration Here ###################

QSMODE = MODE_MPMATH
DPS_MPMATH = 100
DPS_PYTHONTYPES = 25

##########################################################
##########################################################

if QSMODE == MODE_NORM:
    DPS = DPS_PYTHONTYPES
else:
    DPS = DPS_MPMATH
mpmath.mp.dps = DPS

if QSMODE == MODE_NORM:
    QSPI = cmath.pi
else:
    QSPI = mpmath.pi 

############### BASIC TYPES ###############

def QSfloat(val):
    if QSMODE == MODE_NORM:
        return float(val)
    else:
        return mpmath.mpf(val)
    
def QScomplex(val):
    if QSMODE == MODE_NORM:
        return complex(val)
    else:
        if type(val) is str or type(val) is unicode:
            real = None
            imag = None
            delim = None
            if '+' in val[1:]:
                delim = '+'
            elif '-' in val[1:]:
                delim = '-'
            if delim is None:
                if 'j' in val:
                    imag = val.replace('j','')
                else:
                    real = val
            else:
                index = val[1:].find(delim) + 1
                real = val[:index]
                imag = val[index:].replace('j','')
            return mpmath.mpc(real=real,imag=imag)
        else:
            return mpmath.mpc(val)

def QSToSympy(val):
    if QSMODE == MODE_NORM:
        return val
    else:
        return sy.Float(str(val.real),DPS) + sy.Float(str(val.imag),DPS)*sy.I

############### BASIC OPERATIONS ###############

def QSpow(x, y):
        if QSMODE == MODE_NORM:
            return pow(x, y)
        else:
            return mpmath.power(x, y)

def QSexp(x):
        if QSMODE == MODE_NORM:
            return cmath.exp(x)
        else:
            return mpmath.exp(x)

def QSsqrt(x):
        if QSMODE == MODE_NORM:
            return cmath.sqrt(x)
        else:
            return mpmath.sqrt(x)

def QStan(x):
        if QSMODE == MODE_NORM:
            return cmath.tan(x)
        else:
            return mpmath.tan(x)

def QSpolar(x):
        if QSMODE == MODE_NORM:
            return cmath.polar(x)
        else:
            return mpmath.polar(x)
    
############### MATRIX TYPES ###############
    
def QSmatrix(val):
    if QSMODE == MODE_NORM:
        return np.matrix(val, dtype=np.complex128)
    else:
        return mpmath.matrix(val)
    
def QSsqZeros(sz):
    if QSMODE == MODE_NORM:
        return np.matrix(np.zeros((sz, sz), dtype=np.complex128))
    else:
        return mpmath.zeros(sz)
    
def QSidentity(sz):
    if QSMODE == MODE_NORM:
        return np.matrix(np.identity(sz, dtype=np.complex128))
    else:
        return mpmath.eye(sz)

############### MATRIX OPERATIONS ###############
    
def QSshape(mat):
    if QSMODE == MODE_NORM:
        return mat.shape
    else:
        return (mat.rows, mat.cols)

def QSsize(mat):
    if QSMODE == MODE_NORM:
        return mat.size
    else:
        return mat.rows*mat.cols

def QSinvert(mat):
    if QSMODE == MODE_NORM:
        return np.linalg.inv(mat)
    else:
        return mpmath.inverse(mat) 

def QSdot(matA, matB):
    if QSMODE == MODE_NORM:
        return np.dot(matA, matB) 
    else:
        return matA * matB
    
def QSgetRow(mat, m):
    if QSMODE == MODE_NORM:
        return mat[m].tolist()[0]
    else:
        row = []
        for n in range(mat.cols):
            row.append(mat[m,n])
        return row

def QSdet(mat):
    if QSMODE == MODE_NORM:
        return np.linalg.det(mat)
    else:
        return mpmath.det(mat)

def QSsumElements(mat):
    if QSMODE == MODE_NORM:
        XS = 0.0
        for x in np.nditer(mat, flags=['refs_ok']):
            XS += x
    else:
        XS = mpmath.mpc(0.0)
        for i in range(mat.rows):
            for j in range(mat.cols):
                XS += mat[i,j]
    return XS
    
############### OTHER ###############

def formattedFloatString(val, dps):
    if QSMODE == MODE_NORM:
        return ("{:1."+str(dps)+"f}").format(val)
    else:
        return mpmath.nstr(val, mpIntDigits(val)+dps)

def formattedComplexString(val, dps):
    if val.imag < 0.0:
        signStr = ""
    else:
        signStr = "+"
    return formattedFloatString(val.real, dps) + signStr + formattedFloatString(val.imag, dps)+"j"
        
def QSfloatList(lst):
    return str(map(lambda x:str(x),lst)).replace("'","")

def mpIntDigits(num):
    if not mpmath.almosteq(num,0):
        a = mpmath.log(abs(num), b=10)
        b = mpmath.nint(a)
        if mpmath.almosteq(a,b):
            return int(b)+1
        else:
            c = mpmath.ceil(a)
            try:
                return int(c)
            except:
                pass
    else:
        return 0