from __future__ import division
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef outer_and_add_density(np.complex128_t [:] fr1,
                            np.float64_t [:,::1] fthetaphi1,
                            np.complex128_t [:,:,::1] result1,
                            np.complex128_t [:] fr2,
                            np.float64_t [:,::1] fthetaphi2,
                            np.complex128_t [:,:,::1] result2,
                            np.complex128_t [:] fr3,
                            np.float64_t [:,::1] fthetaphi3,
                            np.complex128_t [:,:,::1] result3,
                            np.complex128_t [:] fr4,
                            np.float64_t [:,::1] fthetaphi4,
                            np.complex128_t [:,:,::1] result4,
                            double norm):
    
    cdef unsigned int i,j,k,Nr,Ntheta,Nphi
    cdef double complex fr1i,fr2i,fr3i,fr4i
    
    Nr=fr1.shape[0]
    Ntheta=fthetaphi1.shape[0]
    Nphi=fthetaphi1.shape[1]
    
    for i in xrange(Nr):
        fr1i=fr1[i]
        fr2i=fr2[i]
        fr3i=fr3[i]
        fr4i=fr4[i]
        for j in xrange(Ntheta):
            for k in xrange(Nphi):
                result1[i,j,k]=result1[i,j,k]+fr1i*(fthetaphi1[j,k]*norm)
                result2[i,j,k]=result2[i,j,k]+fr2i*(fthetaphi2[j,k]*norm)
                result3[i,j,k]=result3[i,j,k]+fr3i*(fthetaphi3[j,k]*norm)
                result4[i,j,k]=result4[i,j,k]+fr4i*(fthetaphi4[j,k]*norm)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef outer_and_add_speed(np.complex128_t [:] fr1,
                            np.float64_t [:,::1] fthetaphi1,
                            np.complex128_t [:] fr2,
                            np.float64_t [:,::1] fthetaphi2,
                            np.complex128_t [:,:,::1] result1,
                            np.complex128_t [:] fr3,
                            np.float64_t [:,::1] fthetaphi3,
                            np.complex128_t [:] fr4,
                            np.float64_t [:,::1] fthetaphi4,
                            np.complex128_t [:,:,::1] result2,
                            double norm):
    
    cdef unsigned int i,j,k,Nr,Ntheta,Nphi
    cdef double complex fr1i,fr2i,fr3i,fr4i
    
    Nr=fr1.shape[0]
    Ntheta=fthetaphi1.shape[0]
    Nphi=fthetaphi1.shape[1]
    
    for i in xrange(Nr):
        fr1i=fr1[i]
        fr2i=fr2[i]
        fr3i=fr3[i]
        fr4i=fr4[i]
        for j in xrange(Ntheta):
            for k in xrange(Nphi):
                result1[i,j,k]=result1[i,j,k] +(fr1i*(fthetaphi1[j,k]*norm) +fr2i*(fthetaphi2[j,k]*norm))
                result2[i,j,k]=result2[i,j,k]+(fr3i*(fthetaphi3[j,k]*norm)  +fr4i*(fthetaphi4[j,k]*norm))

