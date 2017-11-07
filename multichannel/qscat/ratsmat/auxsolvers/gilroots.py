# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 20:15:35 2015

@author: gil
@title: Rootfinder

Find the roots of a function f in the complex plane inside of a rectangular region.
We implement the method in the following paper:

 Delves, L. M., and J. N. Lyness. "A numerical method for locating the zeros of
 an analytic function." Mathematics of computation 21.100 (1967): 543-560.

Alternative code using a similar method can be found here:

 http://cpc.cs.qub.ac.uk/summaries/ADKW_v1_0.html

The main idea is to compute contour integrals of functions of the form
:math:`z^k f'/f` around the contour, for integer values of k. Here :math:`f'` denotes the
derivative of f. The resulting values of the contour integrals are proportional
to :math:`\sum_i z_i^k`, where i is the index of the roots.

Throughout we denote :math:`f_{frac} = f'/f`.

I have also tried several optimizations and strategies for numerical stability.

"""
import numpy as np
from scipy import integrate
import math
from gilfunctions import *


mode_default = 0
mode_add_conjs = 0x1

mode_log_recursive = 0x200
mode_log_summary = 0x400
mode_log_all_warn = 0x800
mode_log_debug = 0x1000

#Used for switching the log off on recursion:
mode_log_switch = 0x1FF

warn_imprecise_roots = 1
warn_max_steps_reached = 2
warn_no_bnd_muller_root = 4
warn_bnd_muller_exception = 8
warn_no_int_muller_root = 16
warn_int_muller_exception = 32
warn_not_all_interior_fnd = 64
warn_root_subtraction_division_by_zero = 128

def root_purge(lst,eps=1e-7,min_i=1e-10):
    if len(lst) == 0:
        return []
    for el in lst[:-1]:
        if abs(el-lst[-1]) < eps and \
        (el.imag/lst[-1].imag>=0 or abs(el.imag)<min_i):
            return root_purge(lst[:-1],eps,min_i)
    return root_purge(lst[:-1],eps,min_i) + [lst[-1]]

def add_missing_conjugates(lst,eps=1e-7,min_i=1e-10):
    new_lst = []
    for el in lst:
        new_lst.append(el)
        new_lst.append(el.conjugate())
    return root_purge(new_lst,eps,min_i)

def locate_muller_root(x1,x2,x3,p,roots,failed_roots):
    warn = 0
    try:          
        mull_root,ret = muller(x1,x2,x3,p.f,p.mul_N,p.mul_ltol,p.mul_htol)
        if ret:
            roots.append(mull_root)
        else:
            failed_roots.append(mull_root)
            warn |= p.handle_warning(warn_no_bnd_muller_root)
    except:
        warn |= p.handle_warning(warn_bnd_muller_exception)
    return warn

def correct_roots(roots,rx,ry,rw,rh,add_conjs,dist_eps,min_i):
    roots_inside = inside_boundary(roots,rx,ry,rw,rh)
    conjs_added = 0
    if add_conjs:
        roots_final = inside_boundary(add_missing_conjugates(roots_inside,dist_eps,min_i),
                                      rx,ry,rw,rh)
        conjs_added = len(roots_final)-len(roots_inside)
    else:
        roots_final = roots_inside
    return roots_final, conjs_added

class root_container:
    def __init__(self, known):
        self.known = known
        
        self.boundary_all = []
        self.boundary_failed_mulls = []  
        self.boundary_purged = [] 
        self.boundary_unique = []     
        self.boundary_within = []
        self.boundary_known_unique = []
        self.residues_subtracted = []
        
        self.subtracted = []
        
        self.interior_rough = []
        self.interior_mull_all = []
        self.interior_mull = []
        self.interior_failed_mulls = []
        self.interior_mull_unique = []
        self.interior_mull_final = []
        
        self.interior_all_subs = []
          
        self.all = []
        self.final = []
        self.new = []

    def log(self,warn,roots,conjs_added,I0,lvl_cnt,dist_eps,mode):
        s = " "*lvl_cnt
        d = "-"*lvl_cnt
        if mode&mode_log_summary:
            if warn == 0:
                print d+"Calculations complete. No warnings."
            else:
                print (d+"Calculations completed with following warnings " + 
                        "occurring at least once:")
                if warn & warn_imprecise_roots:
                    print s+"  -Imprecise number of roots in region."
                if warn & warn_max_steps_reached:
                    print s+"  -Number of region steps reached."
                if warn & warn_no_bnd_muller_root:
                    print (s+"  -No boundary muller root found with " + 
                            "specified parameters.")
                if warn & warn_bnd_muller_exception:
                    print s+"  -Exception during Muller routine."
                if warn & warn_no_int_muller_root:
                    print (s+"  -No interior muller root found with " + 
                            "specified parameters.")
                if warn & warn_not_all_interior_fnd:
                    print s+"  -Not all predicted interior roots found."
                if warn & warn_root_subtraction_division_by_zero:
                    print s+"  -Division by zero when subtracting roots."

            num_known_roots = len(self.known)
            num_interior_roots_fnd = len(self.interior_mull_unique)
            num_sub_roots_fnd = len(self.interior_all_subs)
            num_roots_found = num_interior_roots_fnd + \
                              num_sub_roots_fnd + \
                              conjs_added
            if num_known_roots != 0:
                print s + str(num_known_roots) + " known roots."
            print s+"Total of " + str(len(self.new)) + " new roots found."
            print (s+"  " + str(len(self.boundary_within)) + 
                   " from Boundary Muller.")
            print (s+"  Internal: {:.5f}".format(abs(I0)) + 
                   " Roche predicted. " + str(num_roots_found) + " located:")
            print s+"    " + str(num_interior_roots_fnd) + " from Poly Muller."
            print s+"    " + str(num_sub_roots_fnd) + " from subregions."
            if conjs_added != 0:
                print s+"    " + str(conjs_added) + " added conjugates."
                
        if mode&mode_log_debug:
            print "\n"+s+"All:\n" + str(np.array(self.boundary_all))
            print (s+"Unique Purged:\n" + 
                   str(np.array(self.boundary_unique)))
            print s+"Subtracted:\n" + str(np.array(self.subtracted))
            print s+"Failed:\n" + str(np.array(self.boundary_failed_mulls))
            print ""
            print s+"Rough:\n" + str(self.interior_rough)
            print s+"Interior:\n" + str(np.array(self.interior_mull_all))
            print s+"Purged:\n" + str(np.array(self.interior_mull))
            print s+"Unique:\n" + str(np.array(self.interior_mull_unique))
            print s+"Failed:\n" + str(np.array(self.interior_failed_mulls))
            print ""
            print s+"Subs:\n" + str(np.array(self.interior_all_subs))
            print ""
            print s+"All:\n" + str(np.array(self.all))
            print s+"Final:\n" + str(np.array(self.final))
            print s+"New:\n" + str(np.array(self.new))

    def calculate_boundary_roots(self, p, b):
        outliers = find_maxes(map(abs,b.y))
        warn = 0
        for index in outliers:
            warn |= locate_muller_root(b.c[index-2],b.c[index+2],b.c[index]/2,p,
                                       self.boundary_all,
                                       self.boundary_failed_mulls)
        self.boundary_purged = purge(self.boundary_all,p.dist_eps)
        self.boundary_unique = get_unique(self.boundary_purged,self.known,
                                          p.dist_eps)
        self.boundary_within = inside_boundary(self.boundary_unique,
                                               p.rx,p.ry,p.rw,p.rh)
        
        self.boundary_known_unique = purge(self.boundary_unique+self.known,
                                           p.dist_eps)
        # We don't need the roots far outside the boundary
        self.subtracted = inside_boundary(self.boundary_known_unique,
                                          p.rx,p.ry,p.rw+2.,p.rh+2.)
        self.residues_subtracted = residues(b.f_frac,self.subtracted,p.lmt_N,p.lmt_eps)
        
        return warn

class parameters:
    def __init__(self,f,fp,rx,ry,rw,rh,N,outlier_coeff,max_steps,max_order,
           mul_N,mul_ltol,mul_htol,mul_off,dist_eps,lmt_N,lmt_eps,I0_tol,
           mode,min_i,lvl_cnt): 
        self.f = f
        self.fp = fp
        self.rx = rx
        self.ry = ry
        self.rw = rw
        self.rh = rh
        self.N = N
        self.outlier_coeff = outlier_coeff
        self.max_steps = max_steps
        self.max_order = max_order
        self.mul_N = mul_N
        self.mul_ltol = mul_ltol
        self.mul_htol = mul_htol
        self.mul_off = mul_off
        self.dist_eps = dist_eps
        self.lmt_N = lmt_N
        self.lmt_eps = lmt_eps
        self.I0_tol = I0_tol
        self.mode = mode
        self.min_i = min_i
        self.lvl_cnt = lvl_cnt

    def print_region_string(self):
        if self.mode & mode_log_recursive:
            s = "-"*self.lvl_cnt
            print ("\n"+s+"Region(rx,ry,rw,rh): "+str(self.rx)+" "+ \
                   str(self.ry)+" "+str(self.rw)+" "+str(self.rh))
            
    def print_num_regions(self, num_regions):
        if self.mode & mode_log_summary:
            s = "-"*self.lvl_cnt
            print s+"Contains " +str(num_regions) + " subregions:"
        
    def handle_warning(self, warn):
        s = " "*self.lvl_cnt
        imprecise_roots = s+"Warning!! Number of roots may be imprecise for " + \
                          "this N. Increase N for greater precision."
        max_steps_exceeded = s+"Warning!! max_steps exceeded. Some interior " + \
                             "roots might be missing."
        no_bnd_muller_root = s+"Warning!! Boundary Muller failed to converge." 
        bnd_muller_exception = s+"Warning!! Exception during boundary Muller." 
        no_int_muller_root = s+"Warning!! Interior Muller failed to converge." 
        int_muller_exception = s+"Warning!! Exception during interior Muller." 
        not_all_interior_fnd = s+"Warning!! Not all predicted interior roots " + \
                               "found."
        root_subtraction_division_by_zero = s+"Warning!! Division by zero " + \
                                            "during root subtraction."
        
        if self.mode & mode_log_all_warn:
            if warn == warn_imprecise_roots:
                print imprecise_roots
            elif warn == warn_max_steps_reached:
                print max_steps_exceeded
            elif warn == warn_no_bnd_muller_root:
                print no_bnd_muller_root
            elif warn == warn_bnd_muller_exception:
                print bnd_muller_exception
            elif warn == warn_no_int_muller_root:
                print no_int_muller_root
            elif warn == warn_int_muller_exception:
                print int_muller_exception
            elif warn == warn_not_all_interior_fnd:
                print not_all_interior_fnd
            elif warn == warn_root_subtraction_division_by_zero:
                print root_subtraction_division_by_zero
        return warn

class boundary:
    def __init__(self, p):
        self.c = get_boundary(p.rx,p.ry,p.rw,p.rh,p.N)
        self.f_frac = lambda z: p.fp(z)/(2j*np.pi*p.f(z))
        self.y = [self.f_frac(z) for z in self.c]
        self.max_ok = abs(p.outlier_coeff*get_max(self.y))
        

def droots(f,fp,rx,ry,rw,rh,N=10,outlier_coeff=100.,max_steps=5,max_order=10,
           mul_N=400,mul_ltol=1e-12,mul_htol=1e-12,mul_off=1e-5,dist_eps=1e-7,
           lmt_N=10,lmt_eps=1e-3,I0_tol=5e-3,mode=mode_default,min_i=1e-6,
           roots_known=[],lvl_cnt=0):
    '''
    I assume f is analytic with simple (i.e. order one) zeros.

    TODO:
    save values along edges if iterating to a smaller rectangle
    extend to other kinds of functions, e.g. function with non-simple zeros.

    Args:
        f (function): the function for which the roots (i.e. zeros) will be
            found.

        fp (function): the derivative of f.

        rx,ry (floats): The center of the rectangle in the complex
            plane.

        rw,rh (floats): half the width and height of the rectangular
            region.

        N (optional[int]): Number of points to sample per edge

        outlier_coeff (optional[float]): multiplier for coefficient used when 
            subtracting poles to improve numerical stability. 
            See new_f_frac_safe.

        max_step (optional[int]): Number of iterations allowed for algorithm to
            repeat on smaller rectangles.

        mul_tol (optional[float]): muller tolerance.

        mul_N (optional[int]): maximum number of iterations for muller.

        mul_off (optional[float]): muller point offset .

        known roots (optional[list of complex numbers]): Roots of f that are
            already known.

        min_i (optional[boolean]): If function is a polynomial then
            roots will occur in a+ib, a-ib pairs. This options takes this mode
            into account when purging roots that are close to the real axis.
            calculation.

    Returns:
        A list of roots for the function f inside the rectangle determined by
            the values rx,ry,rw and rh. Also a warning and
            number of regions in calculation.
    '''
    p = parameters(f,fp,rx,ry,rw,rh,N,outlier_coeff,max_steps,max_order,mul_N,
                   mul_ltol,mul_htol,mul_off,dist_eps,lmt_N,lmt_eps,I0_tol,
                   mode,min_i,lvl_cnt)
    roots = root_container(roots_known)
    warn = 0
    num_regions = 1
    p.print_region_string()
    
    b = boundary(p)
    warn |= roots.calculate_boundary_roots(p,b)

    y_smooth = []
    for y_el,z_el in zip(b.y,b.c):
        val, ret = new_f_frac_safe(b.f_frac,z_el,roots.residues_subtracted,
                                   roots.subtracted,b.max_ok,y_el,lmt_N,lmt_eps)
        y_smooth.append(val)
        if not ret:
            warn |= p.handle_warning(warn_root_subtraction_division_by_zero)
    I0 = integrate.trapz(y_smooth, b.c)  # Approx number of roots not subtracted
    tot_num_interior_pred = int(math.ceil(abs(I0)-I0_tol))
    
    if I0 < max_order:
        # If there's only a few roots, find them.
        if abs(tot_num_interior_pred-I0)>I0_tol:
            warn |= p.handle_warning(warn_imprecise_roots)
    
        if tot_num_interior_pred == 0:
            roots.final,conjs_added = correct_roots(roots.subtracted,rx,
                                                   ry,rw,rh,mode&mode_add_conjs,
                                                   dist_eps,min_i)
            roots.new = get_unique(roots.final,roots.known,dist_eps)
            roots.log(warn,roots,conjs_added,I0,lvl_cnt,dist_eps,mode)
            return roots.new,warn,num_regions

        roots.interior_rough = locate_poly_roots(y_smooth,b.c,tot_num_interior_pred)
        ##TODO: best way to pick points for Muller method below
        for root in roots.interior_rough:
            warn |= locate_muller_root(root-mul_off,root+mul_off,root,p,
                                       roots.interior_mull_all,
                                       roots.interior_failed_mulls)
        roots.interior_mull = purge(roots.interior_mull_all,dist_eps)
    roots.interior_mull_unique = inside_boundary(get_unique(roots.interior_mull,
                                     roots.boundary_unique,dist_eps),
                                     rx,ry,rw,rh)

    roots.interior_mull_final,conjs_added = correct_roots(roots.interior_mull_unique,
                                                          rx,ry,rw,rh,mode&mode_add_conjs,
                                                          dist_eps,min_i)
    roots.all = purge(roots.boundary_unique+roots.interior_mull_final,dist_eps)
    # Don't count the added conjs at thie stage, just pass them to the subregions.
    # This is because the Roche sometimes does not locate both paired roots.
    all_interior_found = len(roots.interior_mull_unique) >= tot_num_interior_pred
    roche_accurate = abs(tot_num_interior_pred-I0)<I0_tol
    was_subs = False
    # If some interior roots are missed or if there were many roots,
    # subdivide the rectangle and search recursively.
    if (I0>=max_order or not all_interior_found or not roche_accurate) and max_steps!=0:
        was_subs = True
        x_list = [rx - rw / 2.,rx - rw / 2., rx + rw / 2.,rx + rw / 2.]
        y_list = [ry - rh / 2.,ry + rh / 2., ry - rh / 2.,ry + rh / 2.]
        p.print_num_regions(len(x_list))
        new_log = mode if mode&mode_log_recursive else mode&mode_log_switch 
        for x,y in zip(x_list,y_list):
            roots.from_subrectangle,newWarn,new_regions = droots(f,fp,
                x,y,rw/2.,rh/2.,N,outlier_coeff,max_steps-1,max_order,
                mul_N,mul_ltol,mul_htol,mul_off,dist_eps,lmt_N,lmt_eps,I0_tol,
                new_log,min_i,roots.all,lvl_cnt+1)
            warn |= newWarn
            num_regions += new_regions
            roots.interior_all_subs.extend(roots.from_subrectangle)
    elif max_steps == 0:
        warn |= p.handle_warning(warn_max_steps_reached)
    tot_num_interior_found = len(roots.interior_mull_final+roots.interior_all_subs)
    if tot_num_interior_found != tot_num_interior_pred:
        warn |= p.handle_warning(warn_not_all_interior_fnd)
    roots.all = purge(roots.all+roots.interior_all_subs,dist_eps)


    roots.final,conjs_added = correct_roots(roots.all,rx,ry,rw,rh,
                                            mode&mode_add_conjs,
                                            dist_eps,min_i)
    roots.new = get_unique(roots.final,roots.known,dist_eps)
    if was_subs and mode&mode_log_summary: print
    
    roots.log(warn,roots,conjs_added,I0,lvl_cnt,dist_eps,mode)    
    return roots.new,warn,num_regions
