from runargsrange import *
tcp_range_type = argparse.ArgumentParser(description="Two Channel Radial Well - Energy Range File", parents=[tcp_range])
tcp_range_type.add_argument("type_", help="Type to Plot")
args = tcp_range_type.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], eneFactor=ENEFACTOR)
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)

def pstr(value):
    if value < 1e-6:
        return "0.0"
    elif abs(value-1.0) < 1e-6:
        return "1.0"
    else: 
        return str(value)

dEne = (args.eneEnd_-args.eneStart_) / float(args.eneSteps_)
ene = args.eneStart_+args.eneComplex_*1j
    
f = open("res.txt","w")

for i in range(0,args.eneSteps_+1,1):
    eneStr = formattedComplexString(ene, 6)
    try:
        smat.setEnergy(ene)
        if args.type_ == "S":
            f.write(eneStr + "\n" + str(smat)+"\n\n")
        elif args.type_ == "absS":
            f.write(eneStr + "  " + pstr(abs(smat[0][0])) + "    " + pstr(abs(smat[0][1])) + "    " + pstr(abs(smat[1][0])) + "    " + pstr(abs(smat[1][1])) + "\n")
        elif args.type_ == "Unitary":
            if smat.isUnitary():
                f.write(eneStr + "  Unitary\n")
            else:
                f.write(eneStr + "  Not Unitary\n")
        elif args.type_ == "UnitaryOp":
            f.write(eneStr + "\n" + str(smat.uniOp())+"\n\n")
    except DCException as inst:
        f.write(eneStr + str(inst) + "\n")
    ene += dEne

f.close