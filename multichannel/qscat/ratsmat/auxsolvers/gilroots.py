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
from itertools import chain
from scipy import integrate
import math
import cmath as cm

def Muller(x1,x2,x3,f,tol=1e-12,N=400,verbose=False):
    '''
    A method that works well for finding roots locally in the complex plane.
    Uses three points for initial guess, x1,x2,x3.

    Args:
        x1,x2,x3 (complex numbers): initial points for the algorithm.

        f (function): complex valued function for which to find roots.

        tol (optional[float]): tolerance.

        N(optional[int]): maximum number of iterations.

        verbose (optional[boolean]): print warnings.

    Returns:
        estimated root of the function f.
    '''
    n = 0
    x = x3

    if x1 == x2:
        if verbose:
            print "Muller needs x1 and x2 different!!!"
            print "x1 = x2 = ", x1
        return x3

    if x2 == x3:
        if verbose:
            print "Muller needs x2 and x3 different!!!"
            print "x2 = x3 = ", x2
        return x3

    if x1 == x3:
        if verbose:
            print "Muller needs x1 and x3 different!!!"
            print "x1 = x3 = ", x1
        return x3


    while n < N and abs(f(x3))>tol:
        n+=1
        q = (x3 - x2) / (x2 - x1)
        A = q * f(x3) - q*(1.+q)*f(x2)+q**2.*f(x1)
        B = (2.*q+1.)*f(x3)-(1.+q)**2.*f(x2)+q**2.*f(x1)
        C = (1.+q)*f(x3)

        D1 = B+cm.sqrt(B**2-4.*A*C)
        D2 = B-cm.sqrt(B**2-4.*A*C)
        if abs(D1) > abs(D2):
            D = D1
        elif D1 == D2 == 0:
            x = x3
            break
        else: D = D2

        x = x3 - (x3-x2)*2.*C / D

        x1 = x2
        x2 = x3
        x3 = x
        #print x

    converged = True
    if abs(f(x))>=tol:
        converged = False
        if verbose:
            print ("Desired tolerance not reached. Please try different parameters in Muller "
                   "for better results.")
    return converged, x

def residues(f_frac,roots):
    '''
    Finds the resides of :math:`f_{frac} = f'/f` given the location of some roots of f.
    The roots of f are the poles of f_frac.

    Args:
        f_frac (function): a complex.

        roots (a list of complex numbers): the roots of f; poles of f_frac.

    Returns:
        A list of residues of f_frac.

    '''
    return [limit(lambda z: (z-root)*f_frac(z),root) for root in roots]

def limit(f,z0,N=10,eps=1e-3):
    '''
    Takes possibly matrix-valued function f and its simple pole z0 and returns
    limit_{z \to val} f(z). Estimates the value based on N surrounding
    points at a distance eps.

    Args:
        f (function): the function for which the limit will be found.

        z0 (complex number): The value at which the limit is evaluated.

        N (int): number of points used in the estimate.

        eps (optional[float]):
            distance from z0 at which estimating points are placed.

    Returns:
        Limit value (complex):
            The estimated value of :math:`limit_{z -> z_0} f(z)`.

    '''
    t=np.linspace(0.,2.*np.pi*(N-1.)/N,num=N)
    c=np.exp(1j*t)*eps
    try:
        s=sum(f(z0 + c_el) for c_el in c)/float(N)
        return s
    except:
        print "Something went wrong in estimating the limit."
        return

def new_f_frac(f_frac,z0,residues,roots,val=None):
    '''
    Functions that evaluate the f_frac after some roots and their residues are subtracted.
    This function does NOT check to see if there is division by zero of if the
    values become too large.

    We assume here that the poles are of order 1.

    Args:
        f_frac (function): function for which roots will be subtracted.

        z0 (complex number): point where new_f_frac is evaluated.

        residues (list of complex numbers): The corresponding residues to subtract.

        roots (list of complex numbers): The corresponding roots to subtract.

        val (optional[complex number]): We can impose a value f_frac(z0) if we wish.

    Returns:
        The new value of f_frac(z0) once the chosen poles have been subtracted.
    '''
    if val == None:
        val = f_frac(z0)
    for res,root in zip(residues,roots):
        val -= res/(z0-root)
    return val

def new_f_frac_safe(f_frac,z0,residues,roots,max_ok,val=None,verbose=False):
    '''
    Functions that evaluate the f_frac after some roots and their residues are subtracted.
    The safe version checks for large values and division by zero.
    If the value of f_frac(z0) is too large, subtracting the roots of f becomes
    numerically unstable. In this case, we approximate the new function f_frac
    by using the limit function.

    We assume here that the poles are of order 1.

    Args:
        f_frac (function): function for which roots will be subtracted.

        z0 (complex number): point where new_f_frac is evaluated.

        residues (list of complex numbers): The corresponding residues to
            subtract.

        roots (list of complex numbers): The corresponding roots to subtract.

        val (optional[complex number]): We can impose a value f_frac(z0) if
            we wish.

        max_ok (float) Maximum absolute value of f_frac(z0 to use).

        verbose (optional[boolean]): print warnings.

    Returns:
        The new value of f_frac(z0) once the chosen poles have been subtracted.
    '''
    try:
        if val == None:
            val = f_frac(z0)
        if abs(val) < max_ok:
            return new_f_frac(f_frac,z0,residues,roots,val)
        else:
            return limit(lambda z: new_f_frac(f_frac,z,residues,roots),z0)
    except ZeroDivisionError:
        if verbose:
            print 'division by zero in new_f_frac_safe'
        return limit(lambda z: new_f_frac(f_frac,z,residues,roots),z0)

def find_roots(y_smooth,c,num_roots_to_find):
    '''
    given the values y_smooth, locations c, and the number to go up to,
    find the roots using the polynomial trick.

    Args:
        y_smooth (list of complex numbers): points along smoothed-out boundary.
    '''
    p=[0]  ##placeholder
    for i in xrange(1,num_roots_to_find+1):
        p.append(integrate.trapz([el*z**i for el,z in zip(y_smooth,c)],c) )
    e = [1.]
    for k in xrange(1,num_roots_to_find+1):
        s = 0.
        for i in xrange(1,k+1):
            s += (-1.)**(i-1)*e[k-i]*p[i]
        e.append(s / k)
    coeff = [e[k]*(-1.)**(2.*num_roots_to_find-k)
        for k in xrange(0,num_roots_to_find+1)]
    return np.roots(coeff)

def purge(lst,eps=1e-5):
    '''
    Get rid of redundant elements in a list. There is a precision cutoff eps.

    Args:
        lst (list): elements.

        eps (optional[float]): precision cutoff.

    Returns:
        A list without redundant elements.

    '''
    if len(lst) == 0:
        return []
    for el in lst[:-1]:
        if abs(el-lst[-1]) < eps:
            return purge(lst[:-1],eps)
    return purge(lst[:-1],eps) + [lst[-1]]

def linspace(c1,c2,num=50):
    '''
    make a linespace method for complex numbers.

    Args:
        c1,c2 (complex numbers): The two points along which to draw a line.

        num (optional [int]): number of points along the line.

    Returns:
        a list of num points starting at c1 and going to c2.
    '''
    x1 = c1.real
    y1 = c1.imag
    x2 = c2.real*(num-1.)/num+x1*(1.)/num
    y2 = c2.imag*(num-1.)/num+y1*(1.)/num
    return [real+imag*1j for real,imag in zip(np.linspace(x1,x2,num=num),
                                              np.linspace(y1,y2,num=num)) ]


def get_boundary(x_cent,y_cent,width,height,N):
    '''
    Make a rectangle centered at x_cent,y_cent. Find points along this rectangle.
    I use the convention that width/height make up half the dimensions of the rectangle.

    Args:
        x_cent,y_cent (floats): the coordinates of the center of the rectangle.

        width,height (float): The (half) width and height of the rectangle.

        N (int): number of points to use along each edge.

    Returns:
        A list of points along the edge of the rectangle in the complex plane.
    '''
    c1 = x_cent-width+(y_cent-height)*1j
    c2 = x_cent+width+(y_cent-height)*1j
    c3 = x_cent+width+(y_cent+height)*1j
    c4 = x_cent-width+(y_cent+height)*1j
    return  linspace(c1,c2,num=N)+\
            linspace(c2,c3,num=N)+\
            linspace(c3,c4,num=N)+\
            linspace(c4,c1,num=N)


def inside_boundary(roots_near_boundary,x_cent,y_cent,width,height):
    '''
    Takes roots and the specification of a rectangular region
    returns the roots in the interior (and ON the boundary) of the region.

    Args:
        roots_near_boundary (list of complex numbers): roots near the boundary.

        x_cent,y_cent (floats): coordinates of the center of the region.

        width,height (floats): The (half) width of height of the rectangle.

    Returns:
        Roots in the interior and on the boundary of the rectangle.
    '''
    return [root for root in roots_near_boundary if
            x_cent - width <= root.real <= x_cent + width and \
            y_cent - height <= root.imag <= y_cent + height]

def get_max(y):
    '''
    return the :math:`IQR + median` to determine a maximum permissible value to use
    in the numerically safe function new_f_frac_safe.

    '''
    q75, q50, q25 = np.percentile(y, [75 , 50, 25])
    IQR = q75-q25
    return q50+IQR

def find_maxes(y):
    '''
    Given a list of numbers, find the indices where local maxima happen.

    Args:
        y(list of floats).

    Returns:
        list of indices where maxima occur.

    '''
    maxes = []
    for i in xrange(-2,len(y)-2):
        if y[i-1] < y[i] > y[i+1]:
            maxes.append(i)
    return maxes

def root_purge(lst,eps=1e-7,min_imag=1e-10):
    '''
    Get rid of redundant elements in a list. There is a precision cutoff eps.

    Args:
        lst (list): elements.

        eps (optional[float]): precision cutoff.

    Returns:
        A list without redundant elements.

    '''
    if len(lst) == 0:
        return []
    for el in lst[:-1]:
        if abs(el-lst[-1]) < eps and \
        (el.imag/lst[-1].imag>=0 or abs(el.imag)<min_imag):
            return root_purge(lst[:-1],eps,min_imag)
    return root_purge(lst[:-1],eps,min_imag) + [lst[-1]]

def add_conjugates(lst,eps=1e-7,min_imag=1e-10):
    new_lst = []
    for el in lst:
        new_lst.append(el)
        new_lst.append(el.conjugate())
    return root_purge(new_lst,eps,min_imag)

def get_roots_rect_summary(warn,num_final_roots,num_added_conj_roots,roots_near_boundary,
                           I0,num_interior_roots_fnd,num_known_roots,x_cent,y_cent,width,
                           height,num_regions,verbose,summary):
    '''
    Return final roots and optionally prints summary of get_roots_rect.

    Args:
        warn (int): warnings generated during get_roots_rect.
        
        num_final_roots (floats): total roots within region found by the routine.
        
        num_added_conj_roots (floats): number of conjugate roots added.
        
        roots_near_boundary (floats): roots found during the smoothing.
        
        num_known_roots (floats): number of roots previously discovered.
        
        I0 (float): number of roots predicted using Roche.
        
        x_cent,y_cent (floats): coordinates of the center of the region.

        width,height (floats): The (half) width of height of the rectangle.

        num_regions (int): Total number of regions includings outer and subregions.

        verbose (optional[boolean]): print all warnings.

        summary (optional[boolean]): print a summary of warnings at end of
            calculation.
    '''
    if summary:
        print ("For region(x,y,w,h): "+str(x_cent)+" "+str(y_cent)+" "
               +str(width)+" "+str(height)+" (" +str(num_regions-1) + " subregions)")
        if warn == 0:
            print "Calculations complete. No warnings."
        else:
            print "Calculations completed with following warnings occurring at least once:"
            if warn & warn_imprecise_roots:
                print "  -Imprecise number of roots in region."
            if warn & warn_max_steps_exceeded:
                print "  -Number of region steps exceeded."
            if warn & warn_no_muller_root:
                print "  -No muller root found with specified parameters."
            if warn & warn_not_all_interior_fnd:
                print "  -Not all predicted interior roots found."

        roots_within_boundary = inside_boundary(purge(roots_near_boundary),
                                                x_cent,y_cent,width,height)
        print "Total of " + str(num_final_roots) + " roots found."
        if num_known_roots != 0:
            print "  " + str(num_known_roots) + " already discovered."
        print "  " + str(len(roots_within_boundary)) + " during smoothing routine."
        print "  " + str(num_interior_roots_fnd) + " from application of Roche (" + "{:.5f}".format(abs(I0)) + " predicted)."
        if num_added_conj_roots is not None:
            print "  " + str(num_added_conj_roots) + " added conjugates.\n"

def handle_warning(warn, verbose):
    imprecise_roots = "Warning!! Number of roots may be imprecise for " + \
                      "this N. Increase N for greater precision."
    max_steps_exceeded = "Warning!! max_steps exceeded. Some interior " + \
                         "roots might be missing."
    no_muller_root = "Warning!! Muller failed to converge." 
    not_all_interior_fnd = "Warning!! Not all predicted interior roots found."
    
    if verbose:
        if warn == warn_imprecise_roots:
            print imprecise_roots
        elif warn == warn_max_steps_exceeded:
            print max_steps_exceeded
        elif warn == warn_no_muller_root:
            print no_muller_root
        elif warn == warn_not_all_interior_fnd:
            print not_all_interior_fnd
    return warn

warn_imprecise_roots = 1
warn_max_steps_exceeded = 2
warn_no_muller_root = 4
warn_not_all_interior_fnd = 8
def get_roots_rect(f,fp,x_cent,y_cent,width,height,N=10,outlier_coeff=100.,
    max_steps=5,mul_tol=1e-12,mul_N=400,mul_off=1e-5,max_order=10,
    purge_eps=1e-7,conj_min_imag=None,verbose=False,summary=False,
    known_roots=[]):
    '''
    I assume f is analytic with simple (i.e. order one) zeros.

    TODO:
    save values along edges if iterating to a smaller rectangle
    extend to other kinds of functions, e.g. function with non-simple zeros.

    Args:
        f (function): the function for which the roots (i.e. zeros) will be
            found.

        fp (function): the derivative of f.

        x_cent,y_cent (floats): The center of the rectangle in the complex
            plane.

        width,height (floats): half the width and height of the rectangular
            region.

        N (optional[int]): Number of points to sample per edge

        outlier_coeff (optional[float]): multiplier for coefficient used when subtracting
            poles to improve numerical stability. See new_f_frac_safe.

        max_step (optional[int]): Number of iterations allowed for algorithm to
            repeat on smaller rectangles.

        mul_tol (optional[float]): muller tolerance.

        mul_N (optional[int]): maximum number of iterations for muller.

        mul_off (optional[float]): muller point offset .

        known roots (optional[list of complex numbers]): Roots of f that are
            already known.

        conj_min_imag (optional[boolean]): If function is a polynomial then
            roots will occur in a+ib, a-ib pairs. This options takes this mode
            into account when purging roots that are close to the real axis.

        verbose (optional[boolean]): print all warnings and summaries for
            subregions.

        summary (optional[boolean]): print a summary of warnings at end of
            calculation.

    Returns:
        A list of roots for the function f inside the rectangle determined by
            the values x_cent,y_cent,width, and height. Also a warning and
            number of regions in calculation.
    '''

    warn = 0
    num_regions = 1

    c = get_boundary(x_cent,y_cent,width,height,N)
    f_frac = lambda z: fp(z)/(2j*np.pi*f(z))
    y = [f_frac(z) for z in c]

    outliers = find_maxes(map(abs,y))

    roots_near_boundary = []
    cnt = 0
    for outlier_index in outliers:
        try:          
            mull_ret, mull_root = Muller(c[outlier_index-2],c[outlier_index+2],
            (c[outlier_index])/2,f,mul_tol,mul_N,verbose)
            if mull_ret:
                roots_near_boundary.append(mull_root)
            else:
                warn |= handle_warning(warn_no_muller_root,verbose)
        except:
            pass

    subtracted_roots = purge(roots_near_boundary+known_roots)

    ## we don't need the roots far outside the boundary
    subtracted_roots = inside_boundary(subtracted_roots,x_cent,y_cent,
                                       width+2.,height+2.)

    max_ok = abs(outlier_coeff*get_max(y))
    subtracted_residues = residues(f_frac,subtracted_roots)
    y_smooth = [new_f_frac_safe(f_frac,z_el,subtracted_residues,
                                subtracted_roots,max_ok,y_el,verbose)
                                for y_el,z_el in zip(y,c)]
    I0 = integrate.trapz(y_smooth, c)  # approx number of roots not subtracted

    num_roots_interior = int(round(abs(I0)))
    interior_roots = []
    

    ## If there's only a few roots, find them.
    if I0 < max_order:
        if num_roots_interior == 0:
            final_roots = inside_boundary(subtracted_roots,x_cent,y_cent,width,height)
            if conj_min_imag:
                final_roots = inside_boundary(add_conjugates(final_roots),
                                              x_cent,y_cent,width,height)
                num_added_conj_roots = len(final_roots)-len(final_roots)
            else:
                final_roots = final_roots
                num_added_conj_roots = None
            get_roots_rect_summary(warn,len(final_roots),num_added_conj_roots,roots_near_boundary,
                                   num_roots_interior,I0,len(known_roots),x_cent,y_cent,width,
                                   height,num_regions,verbose,summary)
            return final_roots,warn,num_regions
        if abs(num_roots_interior-I0)>0.005:
            warn |= handle_warning(warn_imprecise_roots, verbose)

        rough_roots = find_roots(y_smooth,c,num_roots_interior)

        ##TODO: best way to pick points for Muller method below
        mull_roots = []
        for root in rough_roots:
            mull_ret, mull_root = Muller(root-mul_off,root+mul_off,root,f,
                                         mul_tol,mul_N,verbose)
            if mull_ret:
                mull_roots.append(mull_root)
            else:
                warn |= handle_warning(warn_no_muller_root,verbose)
        interior_roots = purge(mull_roots)
        combined_roots = purge(roots_near_boundary+interior_roots)
    else:
        combined_roots = purge(roots_near_boundary)


    ## if some interior roots are missed or if there were many roots,
    ## subdivide the rectangle and search recursively.
    if I0>=max_order or len(combined_roots) < num_roots_interior and max_steps != 0:
        x_list = [x_cent - width / 2.,x_cent - width / 2.,
                  x_cent + width / 2.,x_cent + width / 2.]
        y_list = [y_cent - height / 2.,y_cent + height / 2.,
                  y_cent - height / 2.,y_cent + height / 2.]
        for x,y in zip(x_list,y_list):           
            roots_from_subrectangle,newWarn,new_regions = get_roots_rect(f,fp,x,y,
                width/2.,height/2.,N,outlier_coeff,max_steps-1,mul_tol,mul_N,mul_off,
                max_order,purge_eps,conj_min_imag,verbose,(verbose and summary),
                known_roots=combined_roots)
            warn |= newWarn
            num_regions += new_regions
            interior_roots = purge(interior_roots+roots_from_subrectangle)
        combined_roots = purge(combined_roots+interior_roots)
    elif max_steps == 0:
        warn |= handle_warning(warn_max_steps_exceeded,verbose)

    if len(interior_roots) != num_roots_interior:
        warn |= handle_warning(warn_not_all_interior_fnd,verbose)

    final_roots = inside_boundary(combined_roots,x_cent,y_cent,width,height)
    if conj_min_imag:
        final_roots_mod = inside_boundary(add_conjugates(final_roots),
                                      x_cent,y_cent,width,height)
        num_added_conj_roots = len(final_roots_mod)-len(final_roots)
    else:
        final_roots_mod = final_roots
        num_added_conj_roots = None
    get_roots_rect_summary(warn,len(final_roots_mod),num_added_conj_roots,roots_near_boundary,
                           I0,len(interior_roots),len(known_roots),x_cent,y_cent,width,
                           height,num_regions,verbose,summary)

    return final_roots_mod,warn,num_regions
