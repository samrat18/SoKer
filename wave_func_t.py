from __future__ import division
import numpy as np
import os
import time
import PLegendre
import wave_form
import matplotlib.pyplot as plt
pi=np.pi
sin=np.sin
cos=np.cos
exp=np.exp
rsun = 6.9598e10
diml = rsun
dimc = 1e6

path='/scratch/samrat/kernel/2015-06-12-l-480-f-1.5to4.7-nnu-1580'
numin=1.5
numax=4.7
ndiv=1580
src_lat=pi/1.5
src_lon=pi/5.5
rcv_lat=pi/2.
rcv_lon=pi/1.5
ellmax=480
nlat=ellmax
nlon=nlat*2

lat=np.atleast_2d(np.linspace(0,pi,nlat)).T
lon=np.atleast_2d(np.linspace(0,pi,nlon))

def derivfd(y):
      dy = np.zeros(np.size(y))
      dy = 0.5 * ( np.roll(y,-1) - np.roll(y,1))
      dy[0] = y[1]-y[0]
      dy[np.size(y)-1] = y[np.size(y)-1]-y[np.size(y)-2]
      return dy


def cosine_chi(lat,lon,lat0,lon0): return (sin(lat0)*sin(lat)*cos(lon-lon0) + cos(lat0)*cos(lat))
def dcos_chi_dtheta(lat,lon,lat0,lon0): return (cos(lat0)*sin(lat) - sin(lat0)*cos(lat)*cos(lon-lon0))
def dcos_chi_dphi(lat,lon,lat0,lon0): return -sin(lat0)*sin(lon-lon0)
LP_atrcv,LP1_atrcv,LP2_atrcv=(np.empty((ellmax+1),dtype=float) for i in range(3))
LP_atrcv,LP1_atrcv,LP2_atrcv=PLegendre.Pl_Pl1_Pl2(ellmax,cosine_chi(rcv_lat,rcv_lon,src_lat,src_lon))
greens_func_rr_atrcv,greens_func_thetar_atrcv,greens_func_phir_atrcv=(np.zeros((ndiv),dtype=float) for i in range(3))
nu = np.linspace( numin, numax, ndiv) * 1e-3
omega=(2*pi)*nu*(diml/dimc)
print 'now computing greens function at receiver'
for omegai in xrange(ndiv):
        t=time.time()
        filename = 'omega-'+str(omegai).zfill(4)+'.npz'
        npzfile = os.path.join(path,filename) 
        gdata=np.load(npzfile)
         
        xi_rr_atrcv=gdata['xisrcatrcv']*LP_atrcv

        p_thetar_atrcv=(gdata['psrcatrcv']*LP1_atrcv)*(dcos_chi_dtheta(rcv_lat,rcv_lon,src_lat,src_lon)/omega[omegai]**2)
        p_phir_atrcv=(gdata['psrcatrcv']*LP1_atrcv)*(dcos_chi_dphi(rcv_lat,rcv_lon,src_lat,src_lon)/omega[omegai]**2)

        greens_func_rr_atrcv[omegai],greens_func_thetar_atrcv[omegai],greens_func_phir_atrcv[omegai]=wave_form.xi_atrcv_f_of_freq(xi_rr_atrcv,
                                                                                                        p_thetar_atrcv,
                                                                                                        p_phir_atrcv
                                                                                                        )
        print 'time taken for omegai',omegai,'is',time.time()-t,'sec'

print 'greens function computed'
print 'now doing fourier transform'
xi_rr_atrcv_func_t=np.fft.fft(greens_func_rr_atrcv)
xi_thetar_atrcv_func_t=np.fft.fft(greens_func_thetar_atrcv)
xi_phir_atrcv_func_t=np.fft.fft(greens_func_phir_atrcv)   

#~ plt.plot(xi_rr_atrcv_func_t)
#~ plt.plot(xi_thetar_atrcv_func_t)
plt.plot(xi_rr_atrcv_func_t)

plt.show()

                                    


