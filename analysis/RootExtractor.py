FILENAMES = ["N=4_S=0_E=1197.dat","N=6_S=0_E=1195.dat","N=8_S=0_E=1197.dat","N=10_S=0_E=1197.dat","N=12_S=0_E=1199.dat","N=14_S=0_E=1196.dat","N=16_S=0_E=1185.dat","N=18_S=0_E=1190.dat","N=20_S=0_E=1197.dat"]
X_CENT = 0.15
Y_CENT = 0.0
WIDTH = 0.14
HEIGHT = 0.3

g = open("stats.dat", 'w')

def write(file, str, tabs=0):
  tabstr = "\t"*tabs + str
  print tabstr
  file.write(tabstr +"\n")

write(g, "stats for region (x,y,w,h): " + str(X_CENT) + ", " + str(Y_CENT) + ", " + str(WIDTH) + ", " + str(HEIGHT))
  
for filename in FILENAMES:
  write(g, filename, 1)

  f = open(filename, 'r')
  
  boundVals = []
  first = True
  for line in f:
    if first or "Complete" in line:
      first = False
      continue
    estr = line[line.rfind(']=')+2:line.rfind('i')]+'j'
    eval = complex(estr)
    if eval.real<X_CENT+WIDTH and eval.real>X_CENT-WIDTH and eval.imag<Y_CENT+HEIGHT and eval.imag>Y_CENT-HEIGHT:
      boundVals.append(eval)
      
  write(g, str(len(boundVals)) +" roots in region", 2)
  nomatchcnt = 0
  rootsstr = ""
  for b1 in boundVals:
    foundmatch = False
    for b2 in boundVals:
      c = b1.conjugate()
      if c.real < b2.real+0.0000001 and c.real > b2.real-0.0000001:
        if c.imag < b2.imag+0.0000001 and c.imag > b2.imag-0.0000001:
          foundmatch = True
          break
    if foundmatch:
      rootsstr += "\t\t\t" + str(b1) + "\n"
    else:
      nomatchcnt += 1
      rootsstr += "\t\t\t###" + b1 + " \n"
  if nomatchcnt==0:
    write(g, "All form a reflected pair.", 2)
  else:
    write(g, str(nomatchcnt) + " without a reflected pair.", 2)
  write(g, rootsstr)
            
  f.close()

g.close()