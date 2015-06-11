from __future__ import division
import numpy as np
import os
import time
from scipy import interpolate
import power_spectrum as ps

coreid=int(os.environ['PBS_VNODENUM'])
pi=np.pi
exp=np.exp
sin=np.sin
cos=np.cos
ndiv=1000
numax = 0.5
numin = 5.5
ellmax=480
nproc=1
rsun = 6.9598e10
diml = rsun
dimc = 1e6
src_lat = pi / 2 
src_lon = pi / 4 
rcv_lat = pi / 2 
rcv_lon = pi / 3

path='/scratch/samrat/kernel/2015-06-10|ell-480-frequency-0.5mHzto5.5mHz-divisions-1000'

cosine_chi_1 = sin(src_lat)*sin(rcv_lat)*cos(rcv_lon-src_lon)
cosine_chi_2 = cos(src_lat)*cos(rcv_lat)
cosine_chi_atrcv = cosine_chi_1 + cosine_chi_2
dcos_chi_dtheta_pow2_atrcv = (cos(src_lat)*sin(rcv_lat) - sin(src_lat)*cos(rcv_lat)*cos(rcv_lon-src_lon))**2.
dcos_chi_dphi_pow2_atrcv= (sin(src_lat)*sin(rcv_lon-src_lon))**2.


def PLegendre(x,lmax):
	Pl=np.zeros(lmax+1)

	if (x==1.0):
		for l in xrange(lmax+1):
			Pl[l]=1.0
	elif (x==-1.0):
		for l in xrange(lmax+1):
			Pl[l]=(-1.0)**l
	else:
		Pl[0]=1.0
		Pl[1]=x
		for  l in xrange(2,lmax+1):
			Pl[l]=((2*l-1)*x*Pl[l-1]-(l-1)*Pl[l-2])/l
	
	return Pl
			
			
nu = np.linspace( numin, numax, ndiv) * 1e-3
omega_per_proc=ndiv//nproc
nu0 = 3.2e-3
sigma = 1e-3

S = (exp(-(nu-nu0)**2. / (2*sigma**2)))
LegendrePsrc_at_rcv=PLegendre(cosine_chi_atrcv,ellmax)

def parallel_on_freq(procid):
	power_spectrum_single_proc_rr=np.zeros((ellmax+1),dtype=float)
	power_spectrum_single_proc_thetar=np.zeros((ellmax+1),dtype=float)
	power_spectrum_single_proc_phir=np.zeros((ellmax+1),dtype=float)
	omega_list_on_single_proc=xrange(procid*ndiv//nproc,(procid+1)*ndiv//nproc)
	for omegai in omega_list_on_single_proc:
		t=time.time()
		filename = 'omega-'+str(omegai).zfill(4)+'.npz'
		npzfile = os.path.join(path,filename) 
		gdata=np.load(npzfile) 
		power_spectrum_over_ell_radial=[]
		power_spectrum_over_ell_angular=[]
		pgrrl,pgthetarl,pgphirl=(np.zeros((ellmax+1),dtype=float) for j in range(3))
		ps.power_spectrum_ell(gdata['xisrcatrcv'].real,
								LegendrePsrc_at_rcv,
								gdata['psrcatrcv'].real,
								pgrrl,
								pgthetarl,
								pgphirl,
								cosine_chi_atrcv)
		#~ Source2=S[omegai]**2.
		power_spectrum_single_proc_rr=np.vstack((power_spectrum_single_proc_rr,pgrrl))
		power_spectrum_single_proc_thetar=np.vstack((power_spectrum_single_proc_thetar,pgthetarl)) 
		power_spectrum_single_proc_phir=np.vstack((power_spectrum_single_proc_phir,pgphirl))
		print 'time to compute Power spectrum at omega  ',omegai,' is',(time.time()-t),'sec'
	return np.delete(power_spectrum_single_proc_rr,0,1),np.delete(power_spectrum_single_proc_thetar,0,1),np.delete(power_spectrum_single_proc_phir,0,1)

power_spectrum_r_r, power_spectrum_theta_r,power_spectrum_phi_r=parallel_on_freq(coreid)
#~ print 'shape of power_spectrum_rr ',power_spectrum_radial_comp.shape, 'and size',power_spectrum_radial_comp.nbytes/1e9,'GB'
#~ print 'shape of power_spectrum_angular ',power_spectrum_angular_comp.shape, 'and size',power_spectrum_angular_comp.nbytes/1e9,'GB'
power_file='procid-{:d}-l-{:}-nu-{:}to{:}'.format(coreid,ellmax,numin,numax)
directory='/scratch/samrat/kernel/power_480_ells/'
os.makedirs(directory)
np.savez(os.path.join(directory,power_file), p_rr=power_spectrum_r_r, p_thetar=power_spectrum_theta_r, p_phir=power_spectrum_phi_r)
