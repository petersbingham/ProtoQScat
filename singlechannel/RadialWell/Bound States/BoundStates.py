import math

#Quick and dirty way of finding solutions for radial well.
#Assumes that we will always approach from the left and won't overshoot.

NUMITERS = 100000

class _Calculator():
  def __init__(self, V, a):
    self.V = V
    self.a = a
  def trig(self, e):
    t = math.tan(e)
    if t!=0:
      return -1.0 * e / t
    else:
      return None
  def circ(self, e):
    return math.sqrt(math.pow(self.rad(),2) - math.pow(e,2))
  def rad(self):
    return self.a*math.sqrt(2.0*self.V)
  def energy(self, e):
    return self.V - math.pow(e,2) / (2*math.pow(self.a,2))
  def k(self, e):
    return math.sqrt(2*(self.energy(e)))
  def K(self, e):
    return math.sqrt(2*(self.V - self.energy(e)))

def _f(num):
  return str(num) + " "
  return "{:8.4f}".format(num)
  
def printStates(V, a, allVals=True):  
  c = _Calculator(V, a)
  
  solutions = []
  de = c.rad() / NUMITERS
  looking = True
  
  for i in range(0,NUMITERS):
    e = i * de
    diff = None
    
    t = c.trig(e)
    if t is not None and t>0:
      diff = c.circ(e) - t
      if diff < 0:
        if looking:
          solutions.append(e-de/2)
          looking = False
      else:
        looking = True
  
  first = _f(V)
    
  if allVals:
    for root in solutions:
      print first+_f(root)+_f(c.energy(root))+_f(c.k(root))+_f(c.K(root))
      first = "        "
  else:
    for root in solutions:
      #first += _f(c.energy(root))
      first += _f(c.k(root))
    print first
  
  