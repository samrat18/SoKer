from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

import pyshtools

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef Pl_Pl1_Pl2(unsigned int lmax,double x):

    '''
    Computes Pl and dPl using shtools' PLegendre_d1
    Computes second derivative using the recursion relation
    d2P(x)=(2l-1)/l *(2dP_{l-1}(x) + x d2P_{l-1}(x) )-(l-1)/l*d2P_{l-2}(x)
    '''
    
    cdef unsigned int l;
    cdef np.ndarray[np.float64_t, ndim=1] Pl,Pl1,Pl2;
    
    Pl,Pl1=pyshtools.PLegendre_d1(lmax,x)
    Pl2=np.zeros(lmax+1)
    
    Pl2[0]=0
    Pl2[1]=0
    
    for l in xrange(2,lmax+1):
        Pl2[l]=((2*l-1)/l * (2*Pl1[l-1]+x*Pl2[l-1])
                -(l-1)/l * Pl2[l-2])
                
    return Pl,Pl1,Pl2


      
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef compute_Pl_and_2_derivs_inplace(unsigned int lmax,np.ndarray[np.float64_t,ndim=2] x_array,
                                        np.ndarray[np.float64_t,ndim=4] LP):
    '''Computes Pl, dPl and d2Pl for each element in x_array, for l=0 to lmax.
    Returns arrays of size lmax+1 rows by x.shape'''
    cdef double x
    cdef unsigned int Nxrows,Nxcols
    cdef Py_ssize_t index
#~     cdef np.ndarray[np.float64_t, ndim=3] Pl,Pl1,Pl2
    
    Nxrows=x_array.shape[0]
    Nxcols=x_array.shape[1]
    
#~     Pl=np.zeros((lmax+1,Nxrows,Nxcols))
#~     Pl1=np.zeros((lmax+1,Nxrows,Nxcols))
#~     Pl2=np.zeros((lmax+1,Nxrows,Nxcols))
    
    for rowindex in xrange(Nxrows):
        for colindex in xrange(Nxcols):
            x=x_array[rowindex,colindex]
#~             Pl[:,rowindex,colindex],Pl1[:,rowindex,colindex],Pl2[:,rowindex,colindex]=Pl_Pl1_Pl2(lmax,x)
            LP[0,:,rowindex,colindex],LP[1,:,rowindex,colindex],LP[2,:,rowindex,colindex]=Pl_Pl1_Pl2(lmax,x)
