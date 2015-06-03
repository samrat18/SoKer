rom __future__ import division
import numpy as np
import PLegendre
import os
import tensorproduct

sin=np.sin
cos=np.cos
pi=np.pi
ndiv = 960
nproc=24*2
corenum = int(os.environ['PBS_VNODENUM'])
src_lat = pi / 2 
src_lon = pi / 2 
rcv_lat = pi / 2 
rcv_lon = pi / 4 

nlat = 480 ;    
nlon = 2 * nlat
ellmax = 500 

 

path = '/scratch/samrat/kernel/greens/'
directory_kSS = '/scratch/samrat/kernel/greens/sound_speed_individual/'

if not os.path.exists(directory_kSS):
    os.makedirs(directory_kSS)

directory_kD = '/scratch/jishnu/kernel/greens/density_individual/'
if not os.path.exists(directory_kD):
    os.makedirs(directory_kD)

nr = np.load(os.path.join(path,'omega-100.npz'))['r'].shape

def lat_lon(nlat,nlon,lat0,lon0):

    lon = np.atleast_2d(np.linspace (0, 2*pi, nlon))
    lat = np.atleast_2d(np.linspace (0, pi,   nlat)).T 
    lat0col = np.ones_like(lat)*lat0

    # Pl(cos_chi) ---- 
    # cos_chi = cos(theta)cos(theta0)+sin(theta)sin(theta0)cos(phi-phi0)
    cosine_chi_1 = sin(lat0)*sin(lat)*cos(lon-lon0)
    cosine_chi_2 = cos(lat0)*cos(lat)
    cosine_chi = cosine_chi_1 + cosine_chi_2
    
    term1 = cos(lat0)*sin(lat)
    term2 = sin(lat0)*cos(lat)*cos(lon-lon0)
    term3 = sin(lat0col)*sin(lon-lon0)
    term4 = sin(lat0)*sin(lat)*cos(lon-lon0)
    term5 = cos(lat0)*cos(lat)
    
    dcos_chi_dtheta = term1 - term2
    dcos_chi_dphi= term3
    
    kss_d2p_dcos = dcos_chi_dtheta**2 + dcos_chi_dphi**2
    kss_dp_dcos = -2*(term4 + term5)
    
    return cosine_chi, dcos_chi_dtheta, dcos_chi_dphi, kss_d2p_dcos, kss_dp_dcos

cosine_chi_src,  dcos_chi_dtheta_src, dcos_chi_dphi_src, kss_d2p_dcos_src, kss_dp_dcos_src = lat_lon( 
                                                                                                    nlat,nlon,src_lat,src_lon)

LP_src_deriv=np.empty((3,ellmax+1,nlat,nlon))
PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_chi_src,LP_src_deriv)

cosine_chi_src =None

cosine_chi_rcv, dcos_chi_dtheta_rcv, dcos_chi_dphi_rcv, kss_d2p_dcos_rcv, kss_dp_dcos_rcv = lat_lon( 
                                                                                                    nlat,nlon,rcv_lat,rcv_lon)

LP_rcv_deriv=np.empty((3,ellmax+1,nlat,nlon))
PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_chi_rcv,LP_rcv_deriv)

cosine_chi_rcv =None

kss_d2p_dcos_src_3d=np.atleast_3d(kss_d2p_dcos_src).transpose(2,0,1)
kss_dp_dcos_src_3d=np.atleast_3d(kss_dp_dcos_src).transpose(2,0,1)

LP_src_deriv[2]*=kss_d2p_dcos_src_3d
LP_src_deriv_kss=LP_src_deriv[1]*kss_dp_dcos_src_3d

kss_d2p_dcos_src_3d=None
kss_dp_dcos_src_3d=None

kss_d2p_dcos_rcv_3d=np.atleast_3d(kss_d2p_dcos_rcv).transpose(2,0,1)
kss_dp_dcos_rcv_3d=np.atleast_3d(kss_dp_dcos_rcv).transpose(2,0,1)

LP_rcv_deriv[2]*=kss_d2p_dcos_rcv_3d
LP_rcv_deriv_kss=LP_rcv_deriv[1]*kss_dp_dcos_rcv_3d

kss_d2p_dcos_rcv_3d=None
kss_dp_dcos_rcv_3d=None


kd_dp_dcos_src_3d = np.atleast_3d(dcos_chi_dtheta_src*dcos_chi_dtheta_rcv
                                     + dcos_chi_dphi_src * dcos_chi_dphi_rcv).transpose(2,0,1) 
dcos_chi_dtheta_src = None
dcos_chi_dtheta_rcv = None
dcos_chi_dphi_src = None
dcos_chi_dphi_rcv = None

LP_src_deriv_kd=LP_src_deriv[1]*kd_dp_dcos_src_3d
LP_rcv_deriv_kd=LP_rcv_deriv[1]

kd_dp_dcos_src_3d = None

def sum_over_l_for_omega(iOmega,lon_indx):
       
    kd_xisrc = np.zeros((nr,nlat),dtype=complex)
    kd_xircv = np.zeros((nr,nlat),dtype=complex)
    kd_psrc = np.zeros((nr,nlat),dtype=complex)
    kd_prcv = np.zeros((nr,nlat),dtype=complex)
    
    kss_src = np.zeros((nr,nlat),dtype=complex)
    kss_rcv = np.zeros((nr,nlat),dtype=complex)
    
    filename = 'omega-{1:d}'.format(iOmega)
    npzfile = os.path.join(path,filename) 
    g_data = np.load(npzfile) 

    for ell in xrange(ellmax):
           
        norm=(2*ell+1)/(4*pi)
        
        tensorproduct.outer_and_add_density(
                        xisrc_denkernel,LP_src_deriv[0,ell,:,lon_indx],kd_xisrc,
                        xircv,LP_rcv_deriv[0,ell,:,lon_indx],kd_xircv,
                        psrc_denkernel,LP_src_deriv_kd[ell,:,lon_indx],kd_psrc,
                        prcv,LP_rcv_deriv_kd[ell,:,lon_indx],kd_prcv,
                        norm
                        )
        tensorproduct.outer_and_add_speed(
                        xisrc_sskernel,LP_src_deriv[0,ell,:,lon_indx],
                        psrc_sskernel,(LP_src_deriv[2,ell,:,lon_indx]+LP_src_deriv_kss[ell,:,lon_indx]),kss_src,
                        xircv_sskernel,LP_rcv_deriv[0,ell,:,lon_indx],
                        prcv_sskernel,(LP_rcv_deriv[2,ell,:,lon_indx]+LP_src_deriv_kss[ell,:,lon_indx]),kss_rcv,
                        norm
                        )
    kd_indv = np.real(kd_xisrc * kd_xircv + kd_psrc * kd_prcv)
    kd_xisrc,kd_xircv,kd_psrc,kd_prcv=None,None,None,None
                 
    kss_indv = np.real(kss_src * kss_rcv)
    kss_src,kss_rcv=None,None
    

    return kkd_indv, kss_indv

def sum_over_omega_for_lon(procid):
    kernel_density_over_lon = []                               
    kernel_sspeed_over_lon = []

    for loni in xrange(procid*nlon//nproc:(procid+1)*nlon//nproc):

        kernel_density_loni = np.zeros((nr,nlat),dtype=complex)
        kernel_sspeed_loni = np.zeros((nr,nlat),dtype=complex)

        for omegai in xrange(ndiv):

            kernel_density_loni_omegai, kernel_sspeed_loni_omgeai = sum_over_l_for_omega(omegai,loni)
    
            kernel_density_loni += kernel_density_loni_omegai
            kernel_speed_loni  += kernel_sspeed_loni_omegai
            kernel_density_loni_omegai,kernel_sspeed_loni_omegai=None,None
        kernel_density_over_lon.append(kernel_density_loni)
        kernel_sspeed_over_lon.append(kernel_sspeed_loni)
        
    return kernel_density_over_lon,kernel_sspeed_over_lon
    
fnameDensity = 'Density-longitude-{:d}to{:d}'.format(procid*nlon//nproc,(procid*nlon//nproc-1),'04')
kDensitylonfile = os.path.join( directory_kD, fnameDensity)

fnameSpeed = 'SoundSpeed-omega-{:d}to{:d}'.format(procid*nlon//nproc,(procid*nlon//nproc-1),'04')
kSpeedlonfile = os.path.join( directory_kSS, fnameSpeed)
kernel_density,kernel_sspeed=sum_over_omega_for_lon(core_num)
np.savez(kDensitynfile,kernel_density=kernel_density)
np.savez(kDensitynfile,kernel_sspeed=kernel_sspeed)
