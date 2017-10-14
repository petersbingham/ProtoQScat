from gilroots import *

def f(x):
  return x*x+4.0*x
def df(x):
  return 2.0*x+4.0

def c(x):
  return x*x*x + 4.0*x + 1.0j*x*x
def dc(x):
  return 3.0*x + 4.0 + 2.0j*x
  
def g(x):
  return 4.0*x*x-7.0*x+(2.0-x)**(0.5)
def dg(x):
  return 8.0*x-7.0-0.5*(2.0-x)**(-0.5)

def h(x):
  return 4.0*x*x+(2.0-x)**(0.5)
def dh(x):
  return 8.0*x-0.5*(2.0-x)**(-0.5)
  
def test_roots(fun, roots):
  for root in roots:
    print fun(root)
    
#print get_roots_rect(f, df, 0, 0, 10, 10)
print get_roots_rect(c, dc, 0, 0, 10, 10)
#print get_roots_rect(g, dg, 1, 0, 1, 1)
#print get_roots_rect(h, dh, 1, 0, 1, 1)