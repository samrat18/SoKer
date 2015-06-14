from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
pi=np.pi
@cython.boundscheck(False)
@cython.wraparound(False)

cpdef xi_atrcv_f_of_freq(np.ndarray[np.float64_t,ndim=1] f1,
                        np.ndarray[np.float64_t,ndim=1] f2,
                        np.ndarray[np.float64_t,ndim=1] f3
                        ):
    cdef unsigned int i,lmax
    cdef double norm,result_1,result_2,result_3
    lmax=f1.shape[0]
    result_1=0.0
    result_2=0.0
    result_3=0.0
    
    for i in xrange (lmax+1):
        norm=(2.0*i+1.0)/(4.0*pi)

        result_1 +=f1[i]*norm
      
        result_2 +=f2[i]*norm
        result_3 +=f3[i]*norm
    return result_1, result_2, result_3
