from __future__ import division
import numpy as np
import os
import time
from scipy import interpolate
import power_spectrum as ps

pi=np.pi
exp=np.exp
sin=np.sin
cos=np.cos

coreid=int(os.environ['PBS_VNODENUM'])
path='/scratch/samrat/kernel/2015-06-12-l-480-f-1.5to4.7-nnu-1580'
ndiv=1580
numin = 1.5
numax = 4.7
nlat=600
ellmax=nlat//2 -1
src_lat = pi / 2 
src_lon = pi / 5.5 
rcv_lat = pi / 2 
rcv_lon = pi / 3.5

nlon=nlat*1
nproc=1

rsun = 6.9598e10
diml = rsun
dimc = 1e6


def PLegendre(x,lmax,nlat,nlon):
        t=time.time()
	Pl=np.zeros((lmax,nlat,nlon))
	Pl1=np.zeros((lmax,nlat,nlon))        
        ps.Pl_Pl1(x,Pl,Pl1)
        print 'legendre_polynomials and derivatives computed in', time.time()-t,'secs'
	return Pl, Pl1


def lat_lon_grid(nlat, nlon, lat0,lon0):
	nlon=2*nlat
	lat=np.atleast_2d(np.linspace(0,pi,nlat)).T
	lon=np.atleast_2d(np.linspace(0,2*pi,nlon))
	lat0col = np.ones_like(lat)*lat0
	cosine_chi_1 = sin(lat0)*sin(lat)*cos(lon-lon0)
	cosine_chi_2 = cos(lat0)*cos(lat)
	cosine_chi = cosine_chi_1 + cosine_chi_2
        term1 = cos(lat0)*sin(lat)
	term2 = sin(lat0)*cos(lat)*cos(lon-lon0)
	term3 = sin(lat0col)*sin(lon-lon0)
	dcos_chi_dtheta =term2-term1
	dcos_chi_dphi= -term3
	return cosine_chi,dcos_chi_dtheta,dcos_chi_dphi

cos_chisrc,dcos_chisrcdtheta,dcos_chisrcdphi=lat_lon_grid(nlat,nlon,src_lat,src_lon)
LegendrePsrc, dLegendrePsrc=PLegendre(cos_chisrc,ellmax+1,nlat,nlon)	

dcos_chisrcdtheta_3d=np.atleast_3d(dcos_chisrcdtheta).transpose(2,0,1)
dcos_chisrcdphi_3d=np.atleast_3d(dcos_chisrcdphi).transpose(2,0,1)
dcos_chisrcdtheta,dcos_chisrcdphi=None , None

print "Multiplying"
G_theta_thetaphi=dcos_chisrcdtheta_3d*dLegendrePsrc
dcos_chisrcdtheta_3d=None
G_phi_thetaphi=dcos_chisrcdphi_3d*dLegendrePsrc
dcos_chisrcdphi_3d,dLegendrePsrc=None,None

nu = np.linspace( numin, numax, ndiv) * 1e-3
omega_per_proc=ndiv//nproc
nu0 = 3.2e-3
sigma = 0.5e-3

cutoff_Source = exp(-(nu-nu0)**2.) / (2*sigma**2)
omega=(2*pi*diml/dimc )*nu


def parallel_on_freq(procid):
	
        print "Inside parallel_on_freq"
        power_spectrum_single_proc_rr=np.zeros((ellmax+1),dtype=float)
	power_spectrum_single_proc_thetar=np.zeros((ellmax+1),dtype=float)
	power_spectrum_single_proc_phir=np.zeros((ellmax+1),dtype=float)
	omega_list_on_single_proc=xrange(procid*ndiv//nproc,(procid+1)*ndiv//nproc)
	
        for omegai in omega_list_on_single_proc:
                		
                t=time.time()
		Source2=cutoff_Source[omegai]**2.
		filename = 'omega-'+str(omegai).zfill(4)+'.npz'
		npzfile = os.path.join(path,filename) 
		gdata=np.load(npzfile) 
		
                psrrl=np.zeros((ellmax+1),dtype=float) 
                gthetarl_sph,gphirl_sph=(np.zeros((nlat,nlon),dtype=float) for i in range(2))
		
                ps.power_spectrum_rr(gdata['xisrcatrcv'][:(ellmax+1)]*Source2,psrrl,
                                        gdata['psrcatrcv'][:(ellmax+1)]*(Source2/(omega[omegai])**2),
                                        G_theta_thetaphi,gthetarl_sph,
                                        G_phi_thetaphi,gphirl_sph,								
                                        )
		
		psthetarl=np.zeros((ellmax+1),dtype=float)
		psphirl=np.zeros((ellmax+1),dtype=float)
		
		ps.power_spectrum_thetaphir(gthetarl_sph,psthetarl,
                                                gphirl_sph,psphirl,
                                                nlat,nlon
                                                )

		power_spectrum_single_proc_rr=np.vstack((power_spectrum_single_proc_rr,psrrl))
		power_spectrum_single_proc_thetar=np.vstack((power_spectrum_single_proc_thetar,psthetarl)) 
		power_spectrum_single_proc_phir=np.vstack((power_spectrum_single_proc_phir,psphirl))
		
                print 'time to compute Power spectrum at omega  ',omegai,' is',(time.time()-t),'sec'
	
        return np.delete(power_spectrum_single_proc_rr,0,1),np.delete(power_spectrum_single_proc_thetar,0,1),np.delete(power_spectrum_single_proc_phir,0,1)

print "Calling function parallel on frequency"
power_spectrum_r_r, power_spectrum_theta_r,power_spectrum_phi_r=parallel_on_freq(coreid)

print 'shape of power_spectrum_rr ', power_spectrum_r_r.shape, 'and size',power_spectrum_r_r.nbytes/1e9,'GB'
print 'shape of power_spectrum_thetar ', power_spectrum_theta_r.shape, 'and size',power_spectrum_theta_r.nbytes/1e9,'GB'
print 'shape of power_spectrum_phir ', power_spectrum_phi_r.shape, 'and size',power_spectrum_theta_r.nbytes/1e9,'GB'

power_file='procid-{:d}-l-{:}-nu-{:}to{:}'.format(coreid,ellmax,numin,numax)
directory='/scratch/samrat/kernel/power_lmax_'+str(ellmax)+'/'
if not os.path.exists(directory): os.makedirs(directory)
np.savez(os.path.join(directory,power_file), p_rr=power_spectrum_r_r, p_thetar=power_spectrum_theta_r, p_phir=power_spectrum_phi_r)