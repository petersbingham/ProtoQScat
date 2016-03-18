import sys
from BoundStates import*

if len(sys.argv) < 5:
  print "Not enough arguments"
else:
  a = float(sys.argv[1])
  start = float(sys.argv[2])
  end = float(sys.argv[3])
  steps = int(sys.argv[4])

V = start
dV = (end-start) / float(steps)

for i in range(0,steps):
  printStates(V, a, False)
  V += dV