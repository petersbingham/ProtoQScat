import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities/General')
import SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.16, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(6,3)

POINTS = 200000

compare_ = 0
if len(sys.argv) < 5:
  print "Not enough arguments"
else:
  fileName_ = sys.argv[1]
  t_ = int(sys.argv[2])
  start_ = complex(sys.argv[3])
  end_ = float(sys.argv[4])
  if len(sys.argv)==6:
    compare_ = int(sys.argv[5])

splitName = fileName_.split("_")
a = splitName[2]
V = splitName[3]
  
with open(fileName_, 'r') as f:
  coeff = {}
  for line in f:
    if line[0] == 'N':
      Nfile = int(line[2:])
      coeff[Nfile] = []
    elif len(line)>1 and line[0]!='*':
      coeff[Nfile].append(complex(line.replace('i','j')))

  
Ss_lastReal = []
Ss_lastImag = []
lastN = None  
      
for N in sorted(coeff.keys()):
  def num(k):
    sum = 0 
    for i in range(len(coeff[N])):
      sum += coeff[N][i] * pow(k,i)
    return sum
      
  def denum(k):
    sum = 0 
    for i in range(len(coeff[N])):
      sum += pow(-1,i) * coeff[N][i] * pow(k,i)
    return sum

  ks = []
  Ss_real = []
  Ss_imag = []  
    
  def calValues(k):
    sval = num(k)/denum(k)
    Ss_real.append(sval.real)
    Ss_imag.append(sval.imag)
    #print str(k) + "  \t" + str(sval)

  def copyValues():
    global Ss_lastReal
    global Ss_lastImag
    global lastN
    Ss_lastReal = [val for val in Ss_real]
    Ss_lastImag = [val for val in Ss_imag]
    lastN = N 
  
  def plotCharts(titleLbl, xLbl):
    name = "Approx S matrix for Radial Well, V0="+'%.1f'%float(V)+", a="+'%.1f'%float(a)+", N="+str(N)+", " + titleLbl
    sp.plotSingle(name,ks,[Ss_real,Ss_imag],xLbl,"",legends=["S.real","S.imag"],path="Results/"+name+".png")
    if lastN:
      if compare_ == 1:
        name = "Approx S matrix for Radial Well, V0="+'%.1f'%float(V)+", a="+'%.1f'%float(a)+", " + titleLbl+", S.real for "+str(lastN)+" & "+str(N)
        sp.plotSingle(name,ks,[Ss_lastReal,Ss_real],xLbl,"",legends=["N="+str(lastN),"N="+str(N)],path="Results/"+name+"_"+str(lastN)+"&"+str(N)+".png")
      elif compare_ == 2:
        name = "Approx S matrix for Radial Well, V0="+'%.1f'%float(V)+", a="+'%.1f'%float(a)+", " + titleLbl+", S.imag for "+str(lastN)+" & "+str(N)
        sp.plotSingle(name,ks,[Ss_lastImag,Ss_imag],xLbl,"",legends=["N="+str(lastN),"N="+str(N)],path="Results/"+name+".png")
 
  if t_ == 2:
    k = start_
    dk = (end_-start_.real) / POINTS
    for i in range(POINTS):
      ks.append(k.real)
      calValues(k)
      k += dk
    plotCharts("k.imag="+str(k.imag), "k.real")
    copyValues()
  elif t_ == 3:
    k = start_
    dk = ((end_-start_.imag) / POINTS) * 1.0j
    for i in range(POINTS):
      ks.append(k.imag)
      calValues(k)
      k += dk
    plotCharts("k.real="+str(k.real), "k.imag")
    copyValues()
