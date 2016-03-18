import cmath

def Kval(k,a,V):
  return cmath.sqrt(k*k + 2.0*V)
  
def num1(k,a,V):
  K = Kval(k,a,V)
  return K + k * cmath.tan(K*a) * 1.0j
  
def num2(k,a,V):
  return cmath.cos(-2.0*k*a) + cmath.sin(-2.0*k*a) * 1.0j
  
def denum(k,a,V):
  K = Kval(k,a,V)
  return K - k * cmath.tan(K*a) * 1.0j

def calValues(k, a, V, Ss_real, Ss_imag):
  sval = (num1(k,a,V)*num2(k,a,V))/denum(k,a,V)
  Ss_real.append(sval.real)
  Ss_imag.append(sval.imag)
  #print str(k) + "  \t" + str(sval)