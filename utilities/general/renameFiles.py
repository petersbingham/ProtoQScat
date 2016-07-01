import os
import os.path
base = u"\\\\?\\E:\\Peter's Documents\\PhD\\Code\\Git\\ProtoQScat\\multichannel\\qscat\\ratsmat\\Results_\\"
for (dirpath, dirnames, filenames) in os.walk(base):
    for idx in range(len(dirnames)):
        if "ROOTS-sympy_det(method berkowitz)" in dirnames[idx]:
          oldname = os.path.join(dirpath, dirnames[idx])
          newname = oldname.replace('),', '),sympy_Poly(),')
          print oldname + "\n" + newname
          os.rename(oldname, newname)
          dirnames[idx] = newname