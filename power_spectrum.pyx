from __future__ import division
import numpy as np
<<<<<<< HEAD
import pyshtools as sh
=======
>>>>>>> 03d9fe0e6ac5c1b70e0c52a7c8f75de2d1ce8ca0
cimport numpy as np
cimport cython
pi=np.pi
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(False)
<<<<<<< HEAD
cpdef Pl_Pl1(np.ndarray[np.float64_t,ndim=2] x,
            np.ndarray[np.float64_t,ndim=3] Pl,
            np.ndarray[np.float64_t,ndim=3] Pl1,
            ):
    cpdef unsigned int i,j,l,nlat,nlon,lmax
    lmax,nlat,nlon=Pl.shape[0],Pl.shape[1],Pl.shape[2]
    
    for i in xrange(nlat):
        
        for j in xrange(nlon):
                                
                if (x[i,j]==1.0) : 
                        
                        for l in xrange(lmax):
                                Pl[l,i,j]=1.0
                        
                elif (x[i,j]==-1.0): 
                        for l in xrange(lmax):
                                Pl[l,i,j]=(-1.0)**l
                        
                else:
                        
                        Pl[0,i,j]=1.0
                        Pl[1,i,j]=x[i,j]
                        for l in xrange(2,lmax):
                                Pl[l,i,j]=((2*l-1)*x[i,j]*Pl[l-1,i,j]-(l-1)*Pl[l-2,i,j])/l
                        
                Pl1[0,i,j]=0.0
                Pl1[1,i,j]=1.0
                for l in xrange(2,lmax):
                        Pl1[l,i,j]=1/l *( (2*l-1.)*( Pl[l-1,i,j]+x[i,j]*Pl1[l-1,i,j] )-(l-1.)*Pl1[l-2,i,j] )

cpdef power_spectrum_rr(np.ndarray[np.float64_t,ndim=1] fr1,
                        np.ndarray[np.float64_t,ndim=1] comp_1,
                        np.ndarray[np.float64_t,ndim=1] fr2,
                        np.ndarray[np.float64_t,ndim=3] flatlonell1,                                        
                        np.ndarray[np.float64_t,ndim=2] comp_2,
                        np.ndarray[np.float64_t,ndim=3] flatlonell2,
                        np.ndarray[np.float64_t,ndim=2] comp_3
                        ):
    cdef unsigned int ellmax, l
    cdef double norm, norm_power

    ellmax=fr1.shape[0]
    
    for l in xrange(ellmax):
        norm_power= 1.0/(2*l+1)
        norm=(2*l+1)/(4*pi)
        comp_1[l]=(fr1[l]**2.)*norm_power
        comp_2 +=(fr2[l]*norm)*flatlonell1[l]
        comp_3 +=(fr2[l]*norm)*flatlonell2[l]
        
cpdef power_spectrum_thetaphir(np.ndarray[np.float64_t,ndim=2] flatlon1,
                                np.ndarray[np.float64_t,ndim=1] out_1,
                                np.ndarray[np.float64_t,ndim=2] flatlon2,
                                np.ndarray[np.float64_t,ndim=1] out_2,
                                unsigned int nlat,
                                unsigned int nlon
                                ):
                                    
    cdef unsigned int ellmax, l, m
    ellmax=out_1.shape[0]
    cdef np.ndarray sphexp1=np.zeros((2,ellmax+1,ellmax+1),dtype=complex)
    cdef np.ndarray sphexp2=np.zeros((2,ellmax+1,ellmax+1),dtype=complex)
    
    
    sphexp1=sh.SHExpandDHC(flatlon1,n=nlat,norm=2,sampling=int(nlon//nlat))
    sphexp2=sh.SHExpandDHC(flatlon2,n=nlat,norm=2,sampling=int(nlon//nlat))
    
    for l in xrange(ellmax):
        
        for m in xrange(l+1):
            
            if m==0:
                out_1[l]+= abs(sphexp1[1,l,m])**2.
                out_2[l]+= abs(sphexp2[1,l,m])**2
            else: 
                out_1[l]+= 2.*abs(sphexp1[1,l,m])**2.
                out_2[l]+= 2.*abs(sphexp2[1,l,m])**2
=======

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
        
        
>>>>>>> 03d9fe0e6ac5c1b70e0c52a7c8f75de2d1ce8ca0
