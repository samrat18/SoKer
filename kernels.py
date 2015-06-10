from __future__ import division
import numpy as np
import PLegendre
import os
import tensorproduct
import time
import errno
import pandas as pd
from time import gmtime, strftime

readme=pd.read_csv('/home/samrat/temp',sep='\t',header=None)
path=readme.iat[0,0]
numin=readme.iat[0,1]
numax=readme.iat[0,2]
ndiv=int(readme.iat[0,3])
ellmax=int(readme.iat[0,4])
nnode=15
ppn=240//nnode
nlat = (ellmax+ppn)-ellmax%ppn    
nlon = 2 * nlat

nodeid = int(os.environ['PBS_NODENUM'])
coreid = int(os.environ['PBS_VNODENUM'])
coreid = (coreid-nodeid*ppn)   #check
omega_per_node=ndiv//nnode
#print 'nodeid and coreid are ',nodeid,coreid
sin=np.sin
cos=np.cos
pi=np.pi
src_lat = pi/2 
src_lon = pi/5.5  
rcv_lat = pi/2 
rcv_lon = pi/3.5 

directory_kSS = '/scratch/samrat/kernel/'+strftime("%Y-%m-%d", gmtime())+'|kss|ell-'+str(ellmax)+'src_lon_pi/5.5_rcv_lon_pi/3.5'+'-frequency-'+str(numin)+'mHzto'+str(numax)+'mHz-divisions-'+str(ndiv)
if not os.path.isdir(directory_kSS):
    try:
        os.makedirs(directory_kSS)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise e
        pass

directory_kD = '/scratch/samrat/kernel/'+strftime("%Y-%m-%d", gmtime())+'|kd|ell-'+str(ellmax)+'src_lon_pi/5.5_rcv_lon_pi/3.5'+'-frequency-'+str(numin)+'mHzto'+str(numax)+'mHz-divisions-'+str(ndiv)
if not os.path.isdir(directory_kD):
    try:
        os.makedirs(directory_kD)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise e
        pass
nr = len(np.load(os.path.join(path,'omega-0100.npz'))['r'])
print nr, nlat, nlon
latfull = np.atleast_2d(np.linspace (0, pi,   nlat)).T
lonfull = np.atleast_2d(np.linspace (0, 2*pi, nlon))

def lat_lon(lat0,lon0,procid):
        
	lat=latfull
	lon=lonfull[:,procid*nlon//ppn:(procid+1)*nlon//ppn]
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

#~ t1=time.time()
cosine_chi_src,  dcos_chi_dtheta_src, dcos_chi_dphi_src, kss_d2p_dcos_src, kss_dp_dcos_src = lat_lon(src_lat,src_lon,coreid)
LP_src_deriv=np.empty((3,ellmax+1,nlat,nlon//ppn))

PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_chi_src,LP_src_deriv)
cosine_chi_src =None

kss_d2p_dcos_src_3d=np.atleast_3d(kss_d2p_dcos_src).transpose(2,0,1)
kss_dp_dcos_src_3d=np.atleast_3d(kss_dp_dcos_src).transpose(2,0,1)
kss_d2p_dcos_src, kss_dp_dcos_src= None,None

LP_src_deriv[2]*=kss_d2p_dcos_src_3d
LP_src_deriv_kss=LP_src_deriv[1]*kss_dp_dcos_src_3d
kss_d2p_dcos_src_3d=None
kss_dp_dcos_src_3d=None

cosine_chi_rcv, dcos_chi_dtheta_rcv, dcos_chi_dphi_rcv, kss_d2p_dcos_rcv, kss_dp_dcos_rcv = lat_lon(rcv_lat,rcv_lon,coreid) 
LP_rcv_deriv=np.empty((3,ellmax+1,nlat,nlon//ppn))

PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_chi_rcv,LP_rcv_deriv)
cosine_chi_rcv =None

kss_d2p_dcos_rcv_3d=np.atleast_3d(kss_d2p_dcos_rcv).transpose(2,0,1)
kss_dp_dcos_rcv_3d=np.atleast_3d(kss_dp_dcos_rcv).transpose(2,0,1)
kss_d2p_dcos_rcv, kss_dp_dcos_rcv=None, None

LP_rcv_deriv[2]*=kss_d2p_dcos_rcv_3d
LP_rcv_deriv_kss=LP_rcv_deriv[1]*kss_dp_dcos_rcv_3d
kss_d2p_dcos_rcv_3d=None
kss_dp_dcos_rcv_3d=None


kd_dp_dcos_src_3d = np.atleast_3d(dcos_chi_dtheta_src*dcos_chi_dtheta_rcv+ dcos_chi_dphi_src * dcos_chi_dphi_rcv).transpose(2,0,1) 
dcos_chi_dtheta_src = None
dcos_chi_dtheta_rcv = None
dcos_chi_dphi_src = None
dcos_chi_dphi_rcv = None

LP_src_deriv_kd=LP_src_deriv[1]*kd_dp_dcos_src_3d
LP_rcv_deriv_kd=LP_rcv_deriv[1]
kd_dp_dcos_src_3d = None

#~ print 'time taken to compute the angular dependence is ',(time.time()-t1)

def sum_over_l_for_omega(iOmega):
	   
	kd_xisrc,kd_xircv,kd_psrc,kd_prcv,kss_src,kss_rcv = (np.zeros((nr,nlat,nlon//ppn),dtype=float) for i in range(6))

	filename = 'omega-'+str(iOmega).zfill(4)+'.npz'
	npzfile = os.path.join(path,filename) 
	gdata=np.load(npzfile)
	gdata=gdata 
	
	#~ t2=time.time() 
	for ell in xrange(ellmax+1):
		#~ t2=time.time()        
		norm=(2.*ell+1.0)/(4*pi)
		
		tensorproduct.outer_and_add_density(
						gdata['xisrc_denkernel'][ell].real,LP_src_deriv[0,ell],kd_xisrc,
						gdata['xircv'][ell].real,LP_rcv_deriv[0,ell],kd_xircv,
						gdata['psrc_denkernel'][ell].real,LP_src_deriv_kd[ell],kd_psrc,
						gdata['prcv'][ell].real,LP_rcv_deriv_kd[ell],kd_prcv,
						norm
						)
		tensorproduct.outer_and_add_speed(
						gdata['xisrc_sskernel'][ell].real,LP_src_deriv[0,ell],
						gdata['psrc_sskernel'][ell].real,(LP_src_deriv[2,ell]+LP_src_deriv_kss[ell]),kss_src,
						gdata['xircv_sskernel'][ell].real,LP_rcv_deriv[0,ell],
						gdata['prcv_sskernel'][ell].real,(LP_rcv_deriv[2,ell]+LP_src_deriv_kss[ell]),kss_rcv,
						norm
						)
		#~ print 'time taken for ell ',ell,' is ',(time.time() - t2),'sec','for iOmega',iOmega,'on proc no. ',coreid
	#~ print 'time taken ',(time.time() - t2),'sec','for iOmega',iOmega,'on proc no. ',coreid	
	kd_indv = kd_xisrc * kd_xircv + kd_psrc * kd_prcv
	kd_xisrc,kd_xircv,kd_psrc,kd_prcv=None,None,None,None
				 
	kss_indv=kss_src * kss_rcv
	kss_src,kss_rcv=None,None

	#~ print kd_indv.shape
	return kd_indv, kss_indv

def sum_over_omega_nodewise_for_truncated_long(node_number,procid):
	kernel_density_lonwise = np.zeros((nr,nlat,nlon//ppn),dtype=float)
	kernel_sspeed_lonwise = np.zeros((nr,nlat,nlon//ppn),dtype=float)
	omega_list_on_present_node=np.arange(node_number,ndiv,ndiv//nnode)
	#omega_list_on_present_node_for_present_proc=np.roll(omega_list_on_present_node,-procid*omega_per_node//ppn)
	for omegai in omega_list_on_present_node:

		kernel_density_lonwise_omegai, kernel_sspeed_lonwise_omegai = sum_over_l_for_omega(omegai)
		
		kernel_density_lonwise += kernel_density_lonwise_omegai
		kernel_sspeed_lonwise += kernel_sspeed_lonwise_omegai
		
		kernel_density_lonwise_omegai,kernel_sspeed_lonwise_omegai=None,None
	return kernel_density_lonwise, kernel_sspeed_lonwise

#~ print kernel_density_lon_lat.shape, kernel_sspeed_lon_lat.shape
#~ print kernel_sspeed_lon_lat.nbytes/1e6 ,'Mb'
#~ print kernel_density_lon_lat.nbytes/1e6 ,'Mb'
#~ return kernel_density_lon,kernel_sspeed_lon
#~ t3=time.time()    
filename_Density = 'Density-nodeid-{:d}-procid-{:d}'.format(nodeid,coreid)
kDensityfile = os.path.join( directory_kD, filename_Density)
filename_SoundSpeed = 'SoundSpeed-nodeid-{:d}-procid-{:d}'.format(nodeid,coreid)
kSoundSpeedfile = os.path.join( directory_kSS, filename_SoundSpeed)

density_kernel,sspeed_kernel=sum_over_omega_nodewise_for_truncated_long(nodeid,coreid)
#~ print 'shapes of density and sound speed kernels are ', density_kernel.shape, sspeed_kernel.shape
#~ print 'sizes of density and sound speed kernels are ', density_kernel.nbytes/1e9,'GB	',sspeed_kernel.nbytes/1e9,'GB'
np.savez(kDensityfile,density_kernel=density_kernel)
np.savez(kSoundSpeedfile,sspeed_kernel=sspeed_kernel)


#~ print 'total time taken for one omega is ',(time.time()-t3),'sec'
