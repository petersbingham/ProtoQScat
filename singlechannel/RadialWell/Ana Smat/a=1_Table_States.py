from AnaSmatFuns import *
import sympy.mpmath as mpm
import cmath

a=1.0

def printAnaRoots(V, fndStates):
  print "\n\nFor V=" + str(V)
  def denumWrap(k):
    return denum(k,a,V)
  for state in fndStates:
    try:
      root = mpm.findroot(denumWrap, state)
      diff = abs(root-state)
      print str(root.real) + "\t" + str(root.imag) + "\t" + str(diff)
    except ValueError:
      print "Not Fnd"
      
printAnaRoots(1.3, [0+0.0674953579513j,
0.00000000000001+0.0664779282583j,
-0.00000000000097+0.0650241884044j,
0.00000000000107+0.0650338242543j,
-0.00000000000124+0.0650338242471j,
-13.8430534529+-2.88606909068j,
-10.6400019225+-2.64286416424j,
-7.39706383908+-2.32104895701j,
-4.04534835279+-1.84131273303j,
4.04534831717+-1.84131277142j,
7.39706377529+-2.32104888766j,
10.6400020094+-2.64286461241j,
13.8431433738+-2.88602645179j])

printAnaRoots(3.4, [-11.0832344543+0.00266984223837j,
-10.4892526403+-2.1815946132j,
-10.48925407+-2.18159496189j,
11.0832344543+0.00266984230236j,
10.4892559068+-2.18159436793j,
10.4892540884+-2.18159486205j,
-13.7272795621+-2.41877748627j,
-10.48925407+-2.18159496189j,
-7.17942445293+-1.87156030651j,
-3.63326414441+-1.42261673949j,
-0.00000001770901+1.45613483867j,
3.63326417088+-1.42261675267j,
7.17942447615+-1.87156033393j,
10.4892540884+-2.18159486205j,
13.7272968816+-2.41879605469j,
16.8993180051+-2.67334326665j])

printAnaRoots(5.5, [-0.00000000000002+1.36069189468j,
0.00000000000002+1.86172325428j,
-0.648701313048+2.01979476093j,
-0.00000004941488+2.33748381979j,
0.63500001055+6.98387861047j,
0.00000000000001+4.75796932778j,
-0.00000000000004+2.63870349514j,
0.648701324208+2.01979477307j,
-0.00000004941488+2.33748381979j,
0.63500001055+6.98387861047j,
-16.8041494433+-2.37027923205j,
-13.5923584505+-2.19107753144j,
13.5923587281+-2.19107603371j,
16.820994237+-2.40099341436j])

printAnaRoots(7.6, [-2.37755480726+-1.1146049474j,
-2.37755481451+-1.11460496251j,
2.37755480862+-1.11460495411j,
2.37755481526+-1.11460496328j,
-16.7013114648+-2.24669112673j,
-13.4499681673+-2.0410789575j,
-10.1240220007+-1.81366820768j,
-6.63472677697+-1.52176406045j,
-0.00001270950444+3.02534693024j,
6.63472680541+-1.52176403851j,
10.1240219791+-1.81366819159j,
13.4499680532+-2.04107884003j,
16.7199063109+-2.25056775976j])

printAnaRoots(9.4, [-1.50673352285+-1.0404047756j,
-1.50673349515+-1.04040478668j,
1.50673352573+-1.04040475892j,
1.50673349811+-1.04040478004j,
-16.6075393489+-2.13615262581j,
-13.3241380302+-1.94410740929j,
-9.95630183476+-1.72016269795j,
-6.37575321183+-1.43442208376j,
0.00003650148399+3.52778615433j,
6.37575318128+-1.43442207151j,
9.95630182275+-1.72016272002j,
13.3241375544+-1.94410717273j,
16.5955823229+-2.16653033752j])

printAnaRoots(9.7, [-1.30461217685+-1.02970378512j,
-1.30461217831+-1.02970380165j,
1.30461216381+-1.02970378873j,
1.30461214624+-1.02970376209j,
-16.5876993247+-2.12577912804j,
-13.3028899768+-1.92989549346j,
-9.92785349406+-1.70649554148j,
-6.33120950767+-1.42171480945j,
0.00100184861303+3.60555216195j,
6.33120955759+-1.42171470216j,
9.92785339119+-1.70649554355j,
13.3028902654+-1.92989523681j])

printAnaRoots(10.0, [-1.06429016432+-1.01939460663j,
-1.06429019259+-1.01939460449j,
1.06429015256+-1.0193946109j,
1.06429014453+-1.01939467283j,
-16.5788337906+-2.10842330571j,
-13.2815676654+-1.91614385223j,
-9.89927058323+-1.69328379394j,
-6.28626557689+-1.40944552621j,
-0.00077412585473+3.6817425476j,
6.28626537335+-1.40944543016j,
9.89927055148+-1.69328384926j,
13.2815680735+-1.91614573287j,
16.5678577014+-2.13106138273j])

printAnaRoots(10.3, [-16.5536138734+-2.09934583623j,
-13.2601762885+-1.90283467937j,
-6.24091740931+-1.39758710174j,
-0.749941057927+-1.00945226714j,
0.00158826761736+3.75554867722j,
0.749941104186+-1.00945232863j,
9.87055507675+-1.6804997131j,
13.2601766899+-1.90283253883j,
16.5542349712+-2.09653247502j])