from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(False)
cpdef genmatrix(double [:] a, double complex[:] b, double complex[:] c, double [:] r):
      cdef unsigned int i,nr
      cdef np.ndarray[np.complex128_t,ndim=2] matrix
      nr =r.shape[0]
      # xuniform[1]-xuniform[0]
      matrix = np.zeros((2* nr, 2*nr),dtype=complex)
      for i in xrange(nr):
            matrix[i,i] = -(a[i] -2/r[i])
            matrix[nr+i,nr+i] = a[i]
            matrix[i,nr + i] = -b[i]
            matrix[nr + i,i] = -c[i]
            if (i==0):
                  matrix[i,i+1] = 1./ (r[i+1]-r[i])
                  matrix[nr+i, nr+i+1] = 1./ (r[i+1]-r[i])
            elif (i>0 and i<7 or i<(nr-1) and i >(nr-8) ):
                  matrix[i,i+1] = 1. / (r[i+1]-r[i-1])
                  matrix[nr+i, nr+i+1] = 1. / (r[i+1]-r[i-1])
                  matrix[i, i-1] = -1./(r[i+1]-r[i-1])
                  matrix[nr+i, nr+i -1] = -1. / (r[i+1]-r[i-1])
            elif (i>6 and i<43 or i<(nr-7)and i>(nr-43)):      
                  matrix[i,i-2] = 1./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[i,i-1] = -8./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[i,i+1] = 8./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[i,i+2] = -1./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])

                  matrix[nr+i,nr+i-2] = 1./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[nr+i,nr+i-1] = -8./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[nr+i,nr+i+1] = 8./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
                  matrix[nr+i,nr+i+2] = -1./(r[i-2]-8*r[i-1]+8*r[i+1]-r[i+2])
            elif (i>42 and i< nr-42):       
                  matrix[i,i-3] = -1. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[i,i-2] =  9. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[i,i-1] = -45./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[i,i+1] = 45./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[i,i+2] = -9./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[i,i+3] = 1. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])

                  matrix[nr+i,nr+i-3] = -1. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[nr+i,nr+i-2] =  9. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[nr+i,nr+i-1] = -45./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[nr+i,nr+i+1] = 45./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[nr+i,nr+i+2] = -9./ (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
                  matrix[nr+i,nr+i+3] = 1. / (-r[i-3]+9*r[i-2]-45*r[i-1]+45*r[i+1]-9*r[i+2]+r[i+3])
            else : 
                  matrix[i,i-1] = -1./(r[i]-r[i-1])
      matrix[nr-1,nr-1] = matrix[nr-1,nr-1]+ 1. / (r[nr-1]-r[nr-2])
      matrix[nr,nr] = matrix[nr,nr]- 1. / (r[1]-r[0])
                  
      return matrix      
