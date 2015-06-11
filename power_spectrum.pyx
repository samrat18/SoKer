from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
pi=np.pi
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(False)

cpdef power_spectrum_ell(np.float64_t [:] fr1,
                    np.float64_t [:] flatlon1,
                    np.float64_t [:] fr2,
                    np.float64_t [:] comp_1,
                    np.float64_t [:] comp_2,
                    np.float64_t [:] comp_3,
                    np.float64_t cosine
                    ):
    cdef unsigned int ellmax, i
    cdef double norm
    cdef float temp
    ellmax=fr1.shape[0]
    
    for i in xrange(ellmax):
        norm= 1.0/(2*i+1)
        comp_1[i]=(fr1[i]**2.)*norm
        if (i<(ellmax-1)):temp=(fr2[i+1]*(2*i+3.0)*(i+1)-cosine*i*(2*i+1.0)*fr2[i])**2. *norm**3.0
        else:temp=(-cosine*i*(2*i+1.0)*fr2[i])**2. *norm**3
        comp_2[i]=temp
        comp_3[i]=temp
        
        
