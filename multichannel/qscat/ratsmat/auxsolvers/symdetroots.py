from auxhelper import *

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
        self._printCalStr("Roots", "calculated", wereLoaded)


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
                if not cmp.complexCompare(abs(root), badRoot):
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