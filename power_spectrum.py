from __future__ import division
import numpy as np
import os
import time
from scipy import interpolate

coreid=int(os.environ['PBS_VNODENUM'])
pi=np.pi
exp=np.exp
sin=np.sin
cos=np.cos
ndiv=960
numax = 4.5
numin = 1.5
ellmax=128
nproc=1
rsun = 6.9598e10
diml = rsun
dimc = 1e6
src_lat = pi / 2 
src_lon = pi / 6 
rcv_lat = pi / 2 
rcv_lon = pi / 4 

path='/scratch/samrat/kernel/greens_7-06-15_10Xdamping/'
codedir="/home/samrat/SoKer"

cosine_chi_1 = sin(src_lat)*sin(rcv_lat)*cos(rcv_lon-src_lon)
cosine_chi_2 = cos(src_lat)*cos(rcv_lat)
cosine_chi_atrcv = cosine_chi_1 + cosine_chi_2

def PLegendre(x,lmax):
	Pl=np.zeros(lmax+1)
	Pl1=np.zeros(lmax+1)
	Pl2=np.zeros(lmax+1)

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
	Pl1[0]=0.0
	Pl1[1]=1.0
	Pl2[0]=0.0
	Pl2[1]=0.0

	for l in xrange(2,lmax+1):
		Pl1[l]=1/l *( (2*l-1.)*( Pl[l-1]+x*Pl1[l-1] ) 
                                                            -(l-1.)*Pl1[l-2] )

		Pl2[l]=1/l *( (2*l-1.)*( 2.*Pl1[l-1]+x*Pl2[l-1])
                                                            -(l-1.)*Pl2[l-2])
	return Pl,Pl1,Pl2
			
			
nu = np.linspace( numin, numax, ndiv) * 1e-3
omega_per_proc=ndiv//nproc
dampingfile=os.path.join(codedir,"m585q.4816")
dmpnu = np.loadtxt(dampingfile,usecols= [2])*1e-6
fwhm = np.loadtxt(dampingfile,usecols= [4])*1e-6 *1e1 #Damping 10x
points = zip(dmpnu, fwhm)
points = sorted(points, key=lambda point: point[0])
dmpnu, fwhm = zip(*points)
damping = interpolate.interp1d(dmpnu, fwhm)
nu = (nu + 1j*damping(nu))
nu0 = 3.2e-3
nu0 = (nu0 + 1j*damping(nu0))
sigma = 1e-3
sigma = (sigma + 1j*damping(sigma))

S = (exp(-(nu-nu0)**2. / (2*sigma**2)))
LegendrePsrc_at_rcv,dLegendrePsrc_at_rcv,d2LegendrePsrc_at_rcv=PLegendre(cosine_chi_atrcv,ellmax)

def parallel_on_freq(procid):
	power_spectrum_single_proc_radial=np.zeros((ellmax+1),dtype=float)
	power_spectrum_single_proc_angular=np.zeros((ellmax+1),dtype=float)
        print power_spectrum_single_proc_radial.shape,power_spectrum_single_proc_angular.shape
	omega_list_on_single_proc=xrange(procid*ndiv//nproc,(procid+1)*ndiv//nproc)
	omgea_list_on_single_proc=np.roll(omega_list_on_single_proc,-procid*omega_per_proc)
	for omegai in omega_list_on_single_proc:
		t=time.time()
		filename = 'omega-'+str(omegai).zfill(4)+'.npz'
		npzfile = os.path.join(path,filename) 
		gdata=np.load(npzfile) 
		power_spectrum_over_ell_radial=[]
                power_spectrum_over_ell_angular=[]
		Source2=np.real(S[omegai])**2.
		for ell in xrange(ellmax+1):
			#~ t_ell=time.time()
			xi_atrcv_dueto_src_sqre=((np.real(gdata['xisrcatrcv'][ell]))*LegendrePsrc_at_rcv[ell])**2. *((2*ell+1)/(4*pi))**2
		        p_atrcv_dueto_src_sqre=((np.real(gdata['psrcatrcv'][ell]))*dLegendrePsrc_at_rcv[ell]**2.)**2. *((2*ell+1)/(4*pi))**2
                        #print 'ell',ell
                        #print 'Pl',LegendrePsrc_at_rcv[ell]
                        #print 'dPl',dLegendrePsrc_at_rcv[ell]
                        #print 'd2Pl',d2LegendrePsrc_at_rcv[ell]  
			power_radial=(xi_atrcv_dueto_src_sqre)*(Source2)
			power_angular=p_atrcv_dueto_src_sqre*Source2
                        power_spectrum_over_ell_radial.append(power_radial)
                        power_spectrum_over_ell_angular.append(power_angular)
			#print 'time to compute P at omega and ell ',omegai, ell,' is',(time.time()-t_ell),'sec'
		power_spectrum_over_ell_radial=np.asarray(power_spectrum_over_ell_radial)
		power_spectrum_over_ell_angular=np.asarray(power_spectrum_over_ell_angular)
                #~ print power_spectrum_over_ell.shape
		power_spectrum_single_proc_radial=np.vstack((power_spectrum_single_proc_radial,power_spectrum_over_ell_radial))
                power_spectrum_single_proc_angular=np.vstack((power_spectrum_single_proc_angular,power_spectrum_over_ell_angular))
		print 'time to compute Power spectrum at omega  ',omegai,' is',(time.time()-t),'sec'
	return np.delete(power_spectrum_single_proc_radial,0,1),np.delete(power_spectrum_single_proc_angular,0,1)

power_spectrum_radial_comp, power_spectrum_angular_comp=parallel_on_freq(coreid)
print 'shape of power_spectrum_radial ',power_spectrum_radial_comp.shape, 'and size',power_spectrum_radial_comp.nbytes/1e9,'GB'
print 'shape of power_spectrum_angular ',power_spectrum_angular_comp.shape, 'and size',power_spectrum_angular_comp.nbytes/1e9,'GB'
power_file='procid-{:d}_p'.format(coreid)
dir='/scratch/samrat/kernel/power/'
os.makedirs(dir)
np.savez(os.path.join(dir,power_file), power_spectrum_radial_comp=power_spectrum_radial_comp, power_spectrum_angular_comp=power_spectrum_angular_comp)

