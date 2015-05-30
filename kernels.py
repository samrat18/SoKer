from __future__ import division
import numpy as np
import PLegendre
import os
from tempfile import mkdtemp
import tensorproduct

tempdir=mkdtemp()

pi=np.pi
ndiv = 960
nproc=24*2
procid = int(os.environ['PBS_VNODENUM'])
src_lat = pi / 2 
src_lon = pi / 2 
rcv_lat = pi / 2 
rcv_lon = pi / 4 

nlat = 1024 ;    
nlon = 2 * nlat ; 
ellmax = 500 ;

 

path = '/scratch/jishnu/kernel/greens'
directory_kSS = '/scratch/samrat/kernel/greens/sound_speed_individual/';

if not os.path.exists(directory_kSS):
    os.makedirs(directory_kSS)
directory_kD = '/scratch/samrat/kernel/greens/density_individual/';

if not os.path.exists(directory_kD):
    os.makedirs(directory_kD)
nr = np.load(os.path.join(path,'omega-100.npz'))['r'].shape

def lat_lon(nlat,nlon,lat0,lon0):

    lon = np.atleast_2d(np.linspace (0, 2*pi, nlon)); 
    lat = np.atleast_2d(np.linspace (0, pi,   nlat)).T; 

    term2 = np.sin(lat0)*np.sin(lat)*np.cos(lon-lon0);

    term1 = np.cos(lat0)*np.cos(lat);

    cosine = term1 + term2;

    u = np.cos(lat0)*np.sin(lat);
    v = np.sin(lat0)*np.cos(lat)*np.cos(lon-lon0);

    c1 = u - v; 
    c11 = c1**2;

    lat0col = np.ones_like(lat)*lat0
    c2= np.sin(lat0col)*np.sin(lon-lon0)
    
    c22 = c2**2;
    
    return lat, lon, cosine, c1, c2, c11, c22

lat, lon, cosine_src, c1_src, c2_src, c11_src, c22_src = lat_lon( nlat,nlon,src_lat,src_lon);

LP_src_deriv=np.empty((3,ellmax+1,nlat,nlon))
LP_src_deriv[:]=PLegendre.compute_Pl_and_2_derivs_2Darray(ellmax,cosine_src)

lat , lon, cosine_rcv, c1_rcv, c2_rcv, c11_rcv, c22_rcv = lat_lon( nlat,nlon,rcv_lat,rcv_lon);

LP_rcv_deriv=np.empty((3,ellmax+1,nlat,nlon))
LP_rcv_deriv[:]=PLegendre.compute_Pl_and_2_derivs_2Darray(ellmax,cosine_rcv)

c1src_plus_c2src_3d=np.atleast_3d(c11_src + c22_src).transpose(2,0,1)

LP_src_deriv[2]*=c1src_plus_c2src_3d

c1src_plus_c2src_3d=None

c1rcv_plus_c2rcv_3d=np.atleast_3d(c11_rcv + c22_rcv).transpose(2,0,1)
LP_rcv_deriv[2]*= c1rcv_plus_c2rcv_3d

c1rcv_plus_c2rcv_3d=None

c_Denlat = c1_src * c1_rcv 
c_Denlon = c2_src * c2_rcv 
c_Denlat_plus_denlon_3d = np.atleast_3d(c_Denlat + c_Denlon).transpose(2,0,1) 

LP_src_deriv[1]*=c_Denlat_plus_denlon_3d

out=np.empty((nr,nlat,nlon),dtype=complex)

def sum_over_l_for_omega(iOmega):
       
    kDensity1 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity2 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity3 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity4 = np.zeros((nr,nlat,nlon),dtype=complex)
    
    kSpeed1 = np.zeros((nr,nlat,nlon),dtype=complex)
    kSpeed2 = np.zeros((nr,nlat,nlon),dtype=complex)
    
    filename = 'omega-{1:d}'.format(iOmega)
    npzfile = os.path.join(path,filename) 
    g_data = np.load(npzfile) 

    for ell in xrange(ellmax):
           
        norm=(2*ell+1)/(4*pi)
        
        kDensity1 += tensorproduct.oneD_and_twoD_to_threeD(
                        g_data['xisrc_denkernel'][ell],LP_src_deriv[0,ell],out)
        
        kDensity2 += tensorproduct.oneD_and_twoD_to_threeD(
                        g_data['xircv'][ell],LP_rcv_deriv[0,ell],out)*norm**2.
        
        kDensity3 += tensorproduct.oneD_and_twoD_to_threeD(
                        g_data['psrc_denkernel'][ell],LP_src_deriv[1,ell],out)
        
        kDensity4 += tensorproduct.oneD_and_twoD_to_threeD(
                        g_data['prcv'][ell], LP_rcv_deriv[1,ell],out)*norm**2.
        
        kSpeed1 += tensorproduct.oneD_and_twoD_to_threeD(
                    g_data['xisrc_sskernel'][ell],LP_src_deriv[0,ell],out)
        kSpeed2 += tensorproduct.oneD_and_twoD_to_threeD(
                    g_data['xircv_sskernel'][ell],LP_rcv_deriv[0,ell],out)
        kSpeed1 += tensorproduct.oneD_and_twoD_to_threeD(
                    g_data['psrc_sskernel'][ell],LP_src_deriv[2,ell],out)
        kSpeed2 += tensorproduct.oneD_and_twoD_to_threeD(
                    g_data['prcv_sskernel'][ell],LP_rcv_deriv[2,ell],out)*norm**2.
         
    kDensityn += (kDensity1 * kDensity2 + kDensity3 * kDensity4);
    kSpeedn +=  (kSpeed1 * kSpeed2) ;
    return kDensityn, kSpeedn

kernel_density = np.zeros((nr,nlat,nlon),dtype=complex)    
kernel_sspeed = np.zeros((nr,nlat,nlon),dtype=complex)

for omegai in xrange(procid*ndiv/nproc : (procid*ndiv/nproc+20)):

    kernel_densityi, kernel_sspeedi = sum_over_l_for_omega(omegai)
    
    kernel_density += kernel_densityi
    kernel_speed  += kernel_densityi

    
fnameDensity = 'Density-omega-{:d}to{:d}'.format(procid+1,procid+21,'03');
kDensitynfile = os.path.join( directory_kD, fnameDensity);

fnameSpeed = 'SoundSpeed-omega-{:d}to{:d}'.format(procid+1,procid+21,'03');
kSpeednfile = os.path.join( directory_kSS, fnameSpeed);

np.savez(kDensitynfile,kernel_density=kernel_density)
np.savez(kDensitynfile,kernel_sspeed=kernel_sspeed)
