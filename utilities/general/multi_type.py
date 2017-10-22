import cmath
import mpmath
import numpy as np
import sympy as sy

MODE_NORM = 0
MODE_MPMATH = 1

##########################################################
################### Configuration Here ###################

MTMODE = MODE_MPMATH
DPS_MPMATH = 100
DPS_PYTHONTYPES = 25

##########################################################
##########################################################

if MTMODE == MODE_NORM:
    DPS = DPS_PYTHONTYPES
else:
    DPS = DPS_MPMATH
mpmath.mp.dps = DPS

if MTMODE == MODE_NORM:
    MTPI = cmath.pi
else:
    MTPI = mpmath.pi 

############### BASIC TYPES ###############

# np.percentile need these overrides.
class MTmpc(mpmath.mpc):
    def __lt__(self, other):
        return self.real < other.real
        
    def __le__(self, other):
        return self.real <= other.real
        
    def __gt__(self, other):
        return self.real > other.real
    
    def __ge__(self, other):
        return self.real >= other.real

def MTfloat(val):
    if MTMODE == MODE_NORM:
        return float(val)
    else:
        return mpmath.mpf(val)
    
def MTcomplex(val):
    if MTMODE == MODE_NORM:
        return complex(val)
    else:
        if type(val) is str or type(val) is unicode:
            if 'nan' in val:
                return mpmath.mpc(real='nan',imag='nan')
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

def MTToSympy(val):
    if MTMODE == MODE_NORM:
        return val
    else:
        return sy.Float(str(val.real),DPS) + sy.Float(str(val.imag),DPS)*sy.I

def MTTompmath(val):
    if MTMODE == MODE_NORM:
        return mpmath.mpc(val.real,val.imag)
    else:
        return mpmath.mpc(real=sy.re(val),imag=sy.im(val))
    
############### BASIC OPERATIONS ###############

def MTlt(val1, val2):
    return val1.real < val2.real
    
def MTle(val1, val2):
    return val1.real <= val2.real

def MTgt(val1, val2):
    return val1.real > val2.real
    
def MTge(val1, val2):
    return val1.real >= val2.real

def MTpercentile(a, q, axis=None, out=None, overwrite_input=False, interpolation='linear', keepdims=False):
    if MTMODE == MODE_NORM:
        return np.percentile(a, q, axis, out, overwrite_input, interpolation, keepdims)
    else:
        return np.percentile(map(lambda v: MTmpc(v), a), q, axis, out, overwrite_input, interpolation, keepdims)

def MTpow(x, y):
    if MTMODE == MODE_NORM:
        return pow(x, y)
    else:
        return mpmath.power(x, y)

def MTexp(x):
    if MTMODE == MODE_NORM:
        return cmath.exp(x)
    else:
        return mpmath.exp(x)

def MTsqrt(x):
    if MTMODE == MODE_NORM:
        return cmath.sqrt(x)
    else:
        return mpmath.sqrt(x)

def MTtan(x):
    if MTMODE == MODE_NORM:
        return cmath.tan(x)
    else:
        return mpmath.tan(x)

def MTpolar(x):
    if MTMODE == MODE_NORM:
        return cmath.polar(x)
    else:
        return mpmath.polar(x)
    
############### MATRIX TYPES ###############
    
def MTmatrix(val):
    if MTMODE == MODE_NORM:
        return np.matrix(val, dtype=np.complex128)
    else:
        return mpmath.matrix(val)
    
def MTsqZeros(sz):
    if MTMODE == MODE_NORM:
        return np.matrix(np.zeros((sz, sz), dtype=np.complex128))
    else:
        return mpmath.zeros(sz)
    
def MTidentity(sz):
    if MTMODE == MODE_NORM:
        return np.matrix(np.identity(sz, dtype=np.complex128))
    else:
        return mpmath.eye(sz)

############# MATRIX CHARACTERISTICS #############
    
def MTshape(mat):
    if MTMODE == MODE_NORM:
        return mat.shape
    else:
        return (mat.rows, mat.cols)

def MTsize(mat):
    if MTMODE == MODE_NORM:
        return mat.size
    else:
        return mat.rows*mat.cols

############### MATRIX OPERATIONS ###############
    
def MTdiagonalise(mat):
    if MTMODE == MODE_NORM:
        w, v = np.linalg.eig(mat)
        P = np.transpose(np.matrix(v, dtype=np.complex128))
        return np.dot(P, np.dot(mat, np.linalg.inv(P)))
    else:
        w, v = mpmath.eig(mat)
        P = mpmath.matrix(v).T
        return P * mat * P**-1

def MTinvert(mat):
    if MTMODE == MODE_NORM:
        return np.linalg.inv(mat)
    else:
        return mpmath.inverse(mat) 

def MTdot(matA, matB):
    if MTMODE == MODE_NORM:
        return np.dot(matA, matB) 
    else:
        return matA * matB
    
def MTgetRow(mat, m):
    if MTMODE == MODE_NORM:
        return mat[m].tolist()[0]
    else:
        row = []
        for n in range(mat.cols):
            row.append(mat[m,n])
        return row

def MTcopyRow(src_mat, dest_mat, m):
    newMat = dest_mat.copy()
    if MTMODE == MODE_NORM:
        for n in range(newMat.shape[1]):
            newMat[m,n] = src_mat[m,n]
    else:
        for n in range(newMat.cols):
            newMat[m,n] = src_mat[m,n]
    return newMat

def MTdet(mat):
    if MTMODE == MODE_NORM:
        return np.linalg.det(mat)
    else:
        return mpmath.det(mat)

def MTsumElements(mat):
    if MTMODE == MODE_NORM:
        XS = 0.0
        for x in np.nditer(mat, flags=['refs_ok']):
            XS += x
    else:
        XS = mpmath.mpc(0.0)
        for i in range(mat.rows):
            for j in range(mat.cols):
                XS += mat[i,j]
    return XS
   
def MTtrace(mat):
    if MTMODE == MODE_NORM:
        return np.trace(mat)
    else:
        t = mpmath.mpc(0.0)
        for i in range(mat.rows):
            t += mat[i,i]
        return t
    
def MTatanElements(mat):
    if MTMODE == MODE_NORM:
        return np.arctan(mat)
    else:
        at = mpmath.matrix(mat.rows, mat.cols)
        for i in range(mat.rows):
            for j in range(mat.cols):
                at[i,j] = mpmath.atan(mat[i,j])
        return at

def _toSymMatrix(mat):
    if MTMODE == MODE_NORM:
        return sy.Matrix(mat)
    else:
        symMat = sy.zeros(mat.rows, mat.cols)
        for r in range(mat.rows):
            for c in range(mat.cols):
                symMat[r,c] = mat[r,c]
        return symMat

def _fromSympyToMpmathMatrix(mat):
    mpMat = mpmath.matrix(mat.shape[0])
    for r in range(mat.shape[0]):
        for c in range(mat.shape[0]):
            mpMat[r,c] = MTTompmath(mat[r,c])
    return mpMat

def MTadjugate(mat):
    symMat = _toSymMatrix(mat)
    return _fromSympyToMpmathMatrix(symMat.adjugate())
    
    
############### OTHER ###############

def formattedFloatString(val, dps):
    if MTMODE == MODE_NORM:
        return ("{:1."+str(dps)+"f}").format(val)
    else:
        return mpmath.nstr(val, mpIntDigits(val)+dps)

def formattedComplexString(val, dps):
    if val.imag < 0.0:
        signStr = ""
    else:
        signStr = "+"
    return formattedFloatString(val.real, dps) + signStr + formattedFloatString(val.imag, dps)+"j"
        
def floatList(lst):
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