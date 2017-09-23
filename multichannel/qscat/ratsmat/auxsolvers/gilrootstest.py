import roots

def f(x):
  return x*x+4*x

def df(x):
  return 2*x+4

def g(x):
  return 4*x*x-7*x+(2-x)**(0.5)

def dg(x):
  return 8*x-7-0.5*(2-x)**(-0.5)

def h(x):
  return 4*x*x+(2-x)**(0.5)

def dh(x):
  return 8*x-0.5*(2-x)**(-0.5)

  
  
def f_roots():
  print Roots.get_roots_rect(f, df, 0, 0, 10, 10)
 
def g_roots():
  roots = Roots.get_roots_rect(g, dg, 1, 0, 1, 1)
  print roots
  for root in roots:
    print g(root)
 
def h_roots():
  roots = Roots.get_roots_rect(h, dh, 1, 0, 1, 1)
  print roots
  for root in roots:
    print h(root)