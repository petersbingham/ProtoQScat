import sys
import os
import time
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')

from runbase import *

from analytical.runargsrange import *
spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Poly Root compare", parents=[tcp_range])
args = spArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, ksigns=[QSfloat(1.0),QSfloat(1.0)], massMult=QSfloat(MASSMULT))
anaSmat = getAnaSmat(args, kCal)
smats = getDiscreteAnaSmats(args)

NUM_FND = 0
TIME = 0

def _findERoot(N, starting):
    global NUM_FND
    global TIME
    try:
        ratSmat = getDecimatedRatSmat(args, smats, kCal, kCal, N, anaSmat, True)
        try:
            startTime = time.clock()
            root = ratSmat.findERoot_Multi(starting)
            print str(N) + "  " + str(root)
            if root is not None:
                NUM_FND += 1
            TIME += time.clock() - startTime
        except ValueError:
            print str(N) + "  None"

    except (DCException, sm.MatException) as inst:
      print str(inst)
      sys.exit()

def _findERootFromActual():
    for i in range(2,33):
        _findERoot(i*2,-0.465713059708835)
_findERootFromActual()
'''     
_findERoot(4,-55.74552346516470+0.00000000000007j)        
_findERoot(6,-1.28809330002954+0.00000000000023j)        
_findERoot(8,-0.46710158556732-0.00000000000004j)        
_findERoot(10,-0.46282974672904+0.00000000000354j)        
_findERoot(12,-0.46571050060371-0.00000000007236j) 
_findERoot(14,-0.46571306416684+0.00000000344248j)        
_findERoot(16,-0.46571305881654+0.00000000026859j)        
_findERoot(18,-0.46571303663291-0.00000001767022j)        
_findERoot(20,-0.46571304111336+0.00000008387897j)   
_findERoot(22,-0.46571302539248+0.00000017195625j)        
_findERoot(24,-0.46571293515146+0.00000031252074j)        
_findERoot(26,-0.46571296833325-0.00000014407311j)        
_findERoot(28,-0.46571285632407-0.00000006209891j)        
_findERoot(30,-0.46570425076196+0.00000300710230j)        
_findERoot(32,-0.46571456185054-0.00000048262638j)        
_findERoot(34,-0.46571712829232+0.00001228786384j)        
_findERoot(36,-0.46571440569712+0.00000179970852j)        
_findERoot(38,-0.46572145780015+0.00000349479808j)        
_findERoot(40,-0.46571834548151+0.00000741456669j)        
_findERoot(42,-0.46571489095852+0.00000386614387j)        
_findERoot(44,-0.46570721719322-0.00003668879879j)        
_findERoot(46,-0.46568587223578-0.00001085776246j)        
_findERoot(48,-0.46570902028484+0.00000122175128j)        
_findERoot(50,-0.46570926561052+0.00000319050473j)        
_findERoot(52,-0.46570008129479-0.00000112900045j)        
_findERoot(54,-0.46571738349316-0.00000109925309j)        
_findERoot(56,-0.46571474911454+0.00000298963017j)        
_findERoot(58,-0.46571363595541+0.00000072195700j)        
_findERoot(60,-0.46572733267133+0.00001917309287j)        
_findERoot(62,-0.46570562176092+0.00000267596220j)        
_findERoot(64,-0.46577735060528-0.00000801125692j)
'''
print NUM_FND  
print TIME