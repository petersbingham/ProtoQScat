from auxhelper import *
import gilroots as pydelves
from globalSettings import *

DELVES_MODE_STANDARD = 0
DELVES_MODE_REFLECT_AXIS = 1
DELVES_MODE_REFLECT_POLE = 2
DELVES_MODE_REFLECT_CHECK = 4
DELVES_MODE_MAINTAIN = 8
DELVES_MODE_RETRY = 16

#v2: In _findRoots function matrix elements with poly(...): matLst[len(matLst)-1].append(poly(val))
class SymDetRoots(AuxHelper):
    def __init__(self, suppressCmdOut, solve_method, sympy_nroots_n, sympy_nroots_maxsteps):
        AuxHelper.__init__(self, suppressCmdOut)
        self.solve_method = solve_method
        self.sympy_detArgs = {'method':'berkowitz'} #'bareis''berkowitz''det_LU'
        self.sympy_polyArgs = {} #'lex''grlex'
        self.numpy_rootsArgs = {}
        self.sympy_nrootsArgs = {'n':sympy_nroots_n, 'maxsteps':sympy_nroots_maxsteps, 'cleanup':True}
        typeStrStart = "v2_sympy_det"+getArgDesc(sy_matrix.det, self.sympy_detArgs)+",sympy_Poly"+getArgDesc(sy_polys.Poly.__new__, self.sympy_polyArgs)
        if self.solve_method == "numpy_roots":
            self.typeStr = typeStrStart+",numpy_roots"+getArgDesc(np.roots, self.numpy_rootsArgs)
        elif self.solve_method == "sympy_nroots":
            self.typeStr = typeStrStart+",sympy_nroots"+getArgDesc(sy_polys.polytools.nroots, self.sympy_nrootsArgs)

    def getRoots(self, mat, k, **args):
        self._startLogAction("getRoots")
        if self.solve_method == "numpy_roots":
            ret = self._getRoots_numpy_roots(mat, k)
        elif self.solve_method == "sympy_nroots":
            ret = self._getRoots_sympy_Poly_nroots(mat, k)
        self._endLogAction("getRoots")
        return ret
      
    def _getRoots_numpy_roots(self, mat, k):
        self._startLogAction("mat.det")
        deter = mat.det(**self.sympy_detArgs)
        self._endLogAction("mat.det")
        
        self._startLogAction("sy_polys.Poly")
        a = sy_polys.Poly(deter, k, **self.sympy_polyArgs)
        self._endLogAction("sy_polys.Poly")
        
        self._startLogAction("all_coeffs")
        coeffs = a.all_coeffs()
        self._endLogAction("all_coeffs")
        
        mappedCoeffs = map(lambda val: complex(val), coeffs)
        
        self._startLogAction("np.roots")
        ret = np.roots(mappedCoeffs, **self.numpy_rootsArgs)
        self._endLogAction("np.roots")
        
        self.printCalStr()
        return ret 
    
    def _getRoots_sympy_Poly_nroots(self, mat, k):
        self._startLogAction("mat.det")
        deter = mat.det(**self.sympy_detArgs)
        self._endLogAction("mat.det")
        
        self._startLogAction("sy_polys.Poly")
        a = sy_polys.Poly(deter, k, **self.sympy_polyArgs) 
        self._endLogAction("sy_polys.Poly") 
        
        self._startLogAction("nroots")
        ret = a.nroots(**self.sympy_nrootsArgs)   
        self._endLogAction("nroots") 
        
        self.printCalStr()
        return ret 
    
    def printCalStr(self, wereLoaded=False):
        self._printCalStr("SymDetRoots", "calculated", wereLoaded)


class DelvesRoots(AuxHelper):
    def __init__(self, suppressCmdOut, rx, ry, rw, rh, N, outlier_coeff, max_steps, max_order, 
                 mul_N, mul_ltol, mul_htol, mul_off, dist_eps, lmt_N, lmt_eps, min_i, log, mode):
        AuxHelper.__init__(self, suppressCmdOut)
        self.mode = mode
        self.delves_Args = {'rx':rx, 'ry':ry, 'rw':rw, 'rh':rh, 'N':N, 'outlier_coeff':outlier_coeff, 'max_steps':max_steps, 
                            'mul_N':mul_N, 'mul_ltol':mul_ltol, 'mul_htol':mul_htol, 'mul_off':mul_off, 'max_order':max_order, 
                            'dist_eps':dist_eps, 'lmt_N':lmt_N, 'lmt_eps':lmt_eps, 'min_i':min_i, 'log':log}
        self.typeStr = "droots"+getArgDesc(pydelves.droots, self.delves_Args, ["roots_known", "lvl_cnt"]) + ", mode " +str(mode)
        self.cmp = num.Compare(dist_eps)

    def _doesRootAlreadyExist(self, roots, newRoot):
        found = False
        for root in roots:
            if self.cmp.complexCompare(newRoot, root): 
                found = True
                break
        return found
    
    def _doesRootExist(self, root, f):
        offset = self.delves_Args['muller_offset']
        checkRoot = Muller(root-offset, root+offset, root, f)
        if self.cmp.complexCompare(root, checkRoot):
            return checkRoot
        return None
    
    def _handleNewRoot(self, newRoots, newRoot, f):
        if self.mode & DELVES_MODE_REFLECT_CHECK:
            checkRoot = self._doesRootExist(newRoot, f)
            if checkRoot is not None:
                newRoots.append(checkRoot)
        else:
            newRoots.append(newRoot)
                                                
    def getRoots(self, f, fp, lastRoots):
        if self.mode & DELVES_MODE_REFLECT_AXIS:
            min = self.delves_Args['y_cent'] - self.delves_Args['height']
            max = self.delves_Args['y_cent'] + self.delves_Args['height']
            if abs(max) >= abs(min):
                self.delves_Args['y_cent'] = max / 2.0
                self.delves_Args['height'] = max / 2.0
            else:
                self.delves_Args['y_cent'] = -min / 2.0
                self.delves_Args['height'] = -min / 2.0

        warn = 0x100
        while warn!=0 and (warn==0x100 or self.mode & DELVES_MODE_RETRY):
            if warn & pydelves.warn_imprecise_roots:
                self.delves_Args['N'] *= 2
            if warn & pydelves.warn_max_steps_reached:
                self.delves_Args['max_steps'] *= 2
            if warn & pydelves.warn_no_bnd_muller_root or warn & pydelves.warn_no_int_muller_root:
                self.delves_Args['mul_N'] *= 2
            roots, warn, num_regions = pydelves.droots(f, fp, **self.delves_Args)

        if self.mode & DELVES_MODE_REFLECT_AXIS:
            newRoots = []
            for root in roots:
                newRoot = root.conjugate()
                if (abs(max) >= abs(min) and abs(newRoot.imag)<=abs(min)) or (abs(max) < abs(min) and abs(newRoot.imag)<=abs(max)):
                    self._handleNewRoot(newRoots, newRoot, f)
            roots.extend(newRoots)

        if self.mode & DELVES_MODE_REFLECT_POLE:
            newRoots = []
            for root in roots:
                newRoot = root.conjugate()
                if not self._doesRootAlreadyExist(roots, newRoot):
                    self._handleNewRoot(newRoots, newRoot, f)
            roots.extend(newRoots)           

        if self.mode & DELVES_MODE_MAINTAIN:
            for lastRoot in lastRoots:
                if not self._doesRootAlreadyExist(roots, lastRoot):
                    checkRoot = Muller(lastRoot-offset, lastRoot+offset, lastRoot, f)
                    if self.cmp.complexCompare(checkRoot, lastRoot):
                        roots.append(checkRoot)

        return roots

    def printCalStr(self, wereLoaded=False):
        self._printCalStr("DelvesRoots", "calculated", wereLoaded)
    

class RootClean(AuxHelper):
    def __init__(self, suppressCmdOut, clean_width):
        AuxHelper.__init__(self, suppressCmdOut)
        self.clean_width = clean_width
        if self.clean_width is not None:
            self.typeStr = "dk" + str(self.clean_width) + "," + str(self.clean_width) + "i"
        else:
            self.typeStr = None

    def cleanRoots(self, roots, badRoot):
        self._startLogAction("Cleaning roots")
        cleanRoots = []
        rejectRoots = []
        if self.typeStr is not None:
            cmp = num.Compare(self.clean_width) 
            for root in roots:
                if not cmp.complexCompare(abs(root[k_INDEX]), badRoot):
                    cleanRoots.append(root)
                else:
                    rejectRoots.append(root)
        else:
            cleanRoots = roots
        self.printCalStr()
        self._endLogAction("Cleaning roots") 
        return cleanRoots, rejectRoots
        
    def printCalStr(self, wereLoaded=False):
        self._printCalStr("Roots", "cleaned", wereLoaded)