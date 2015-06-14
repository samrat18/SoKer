from __future__ import division
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef outer_and_add_density(np.float64_t [:] fr1,
                            np.float64_t [:,::1] fthetaphi1,
                            np.float64_t [:,:,::1] result1,
                            np.float64_t [:] fr2,
                            np.float64_t [:,::1] fthetaphi2,
                            np.float64_t [:,:,::1] result2,
                            np.float64_t [:] fr3,
                            np.float64_t [:,::1] fthetaphi3,
                            np.float64_t [:,:,::1] result3,
                            np.float64_t [:] fr4,
                            np.float64_t [:,::1] fthetaphi4,
                            np.float64_t [:,:,::1] result4,
                            np.float64_t [:] fr5,                      
                            np.float64_t [:] fr6,
                            np.float64_t [:,::1] fthetaphi6,
                            np.float64_t [:,:,::1] result5,
                            np.float64_t [:] fr7,
                            np.float64_t [:] fr8,
                            np.float64_t [:,::1] fthetaphi8,
                            np.float64_t [:,:,::1] result6,
                            np.float64_t [:] fr9,
                            np.float64_t [:,::1] fthetaphi9,
                            np.float64_t [:,:,::1] result7,
                            np.float64_t [:] fr10,
                            np.float64_t [:] fr11,
                            np.float64_t [:,::1] fthetaphi11,
                            np.float64_t [:,:,::1] result8,
                            double norm
                            ):
    cdef unsigned int i,j,k, Nr,Ntheta,Nphi
    cdef double fr1i,fr2i,fr3i,fr4i,fr5i,fr6i,fr7i,fr8i
    
    Nr=fr1.shape[0]
    Ntheta=fthetaphi1.shape[0]
    Nphi=fthetaphi1.shape[1]
    
    for i in xrange(Nr):
        fr1i=fr1[i]
        fr2i=fr2[i]
        fr3i=fr3[i]
        fr4i=fr4[i]
        
        fr5i=fr5[i]
        fr6i=fr6[i]
        fr7i=fr7[i]
        fr8i=fr8[i]
        
        fr9i=fr9[i]
        fr10i=fr10[i]
        fr11i=fr11[i]
        
        for j in xrange(Ntheta):
            for k in xrange(Nphi):
                
                #term1
                result1[i,j,k]=result1[i,j,k]+(fr1i*norm)*fthetaphi1[j,k]
                result2[i,j,k]=result2[i,j,k]+(fr2i*norm)*fthetaphi2[j,k]
                result3[i,j,k]=result3[i,j,k]+(fr3i*norm)*fthetaphi3[j,k]
                result4[i,j,k]=result4[i,j,k]+(fr4i*norm)*fthetaphi4[j,k]
                
                #term2
                result5[i,j,k]=result5[i,j,k] +((fr5i*norm)*fthetaphi1[j,k]  +(fr6i*norm)*fthetaphi6[j,k])
                result6[i,j,k]=result6[i,j,k] +((fr7i*norm)*fthetaphi2[j,k]  +(fr8i*norm)*fthetaphi8[j,k])
                
                #term3
                result7[i,j,k]=result7[i,j,k] +((fr9i*norm)*fthetaphi9[j,k])
                result8[i,j,k]=result8[i,j,k] +((fr10i*norm)*fthetaphi2[j,k]  +(fr11i*norm)*fthetaphi11[j,k])
                


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef outer_and_add_speed(np.float64_t [:] fr1,
                            np.float64_t [:,::1] fthetaphi1,
                            np.float64_t [:] fr2,
                            np.float64_t [:,::1] fthetaphi2,
                            np.float64_t [:,:,::1] result1,
                            np.float64_t [:] fr3,
                            np.float64_t [:,::1] fthetaphi3,
                            np.float64_t [:] fr4,
                            np.float64_t [:,::1] fthetaphi4,
                            np.float64_t [:,:,::1] result2,
                            double norm):
    
    cdef unsigned int i,j,k,Nr,Ntheta,Nphi
    cdef double fr1i,fr2i,fr3i,fr4i
    
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
                result1[i,j,k]=result1[i,j,k] +((fr1i*norm)*fthetaphi1[j,k]  +(fr2i*norm)*fthetaphi2[j,k])
                result2[i,j,k]=result2[i,j,k] +((fr3i*norm)*fthetaphi3[j,k]  +(fr4i*norm)*fthetaphi4[j,k])
