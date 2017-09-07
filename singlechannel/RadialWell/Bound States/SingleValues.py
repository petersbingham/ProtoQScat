import sys
from BoundStates import*

if len(sys.argv) < 3:
  print "Not enough arguments"
else:
  V = float(sys.argv[1])
  a = float(sys.argv[2])

printStates(V, a)