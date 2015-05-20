from __future__ import division
import numpy as np
import PLegendre
from scipy.io import loadmat
import os,time,multiprocessing,sys
from tempfile import mkdtemp
import tensorproduct
import thread

node_num=0
#~ node_num = int(os.environ['PBS_NODENUM'])

tempdir=mkdtemp()

pi=np.pi

src_lat = pi / 2 
src_lon = pi / 2 
rcv_lat = pi / 2 
rcv_lon = pi / 2 

nlat = 300 ;    
nlon = 2 * nlat ; 
ellmax = 250 ;
ellmin = 0 ;
nomegadown = 1 ; 
nomegaup = 4 ;
nr = 243 ; 

path = '/home/samrat/final_codes/data_xi_p_21-01-14'
directory_kSS = '/home/samrat/final_codes/12_05_15_kernels_python/sound_speed_individual/';

if not os.path.exists(directory_kSS):
    os.makedirs(directory_kSS)
    
directory_kD = '/home/samrat/final_codes/12_05_15_kernels_python/density_individual/';
if not os.path.exists(directory_kD):
    os.makedirs(directory_kD)

def lat_lon(nlat,nlon,lat0,lon0):

    lon,lat = np.atleast_2d(np.linspace(0, 2*pi, nlon,endpoint=False),np.linspace(0, pi,nlat))
    lat=lat.T

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
    
    return cosine, c1, c2, c11, c22


#~ Source 
cosine_src, c1_src, c2_src, c11_src, c22_src = lat_lon( nlat,nlon,src_lat,src_lon);

LP_src_deriv=np.empty((3,ellmax+1,nlat,nlon))
PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_src,LP_src_deriv)

cosine_src=None

c1src_plus_c2src_3d=np.atleast_3d(c11_src + c22_src).transpose(2,0,1)
c11_src=None
c22_src=None
LP_src_deriv[2]*=c1src_plus_c2src_3d
c1src_plus_c2src_3d=None



#~ Receiver
cosine_rcv, c1_rcv, c2_rcv, c11_rcv, c22_rcv = lat_lon( nlat,nlon,rcv_lat,rcv_lon);

LP_rcv_deriv=np.empty((3,ellmax+1,nlat,nlon))
PLegendre.compute_Pl_and_2_derivs_inplace(ellmax,cosine_rcv,LP_rcv_deriv)

cosine_rcv=None

c1rcv_plus_c2rcv_3d=np.atleast_3d(c11_rcv + c22_rcv).transpose(2,0,1)
c11_rcv=None
c22_rcv=None
LP_rcv_deriv[2]*= c1rcv_plus_c2rcv_3d
c1rcv_plus_c2rcv_3d=None


c_Denlat = c1_src * c1_rcv ;
c_Denlon = c2_src * c2_rcv ;
c1_src=None
c1_rcv=None
c2_src=None
c2_rcv=None
c_Denlat_plus_denlon_3d = np.atleast_3d(c_Denlat + c_Denlon).transpose(2,0,1) ; 

LP_src_deriv[1]*=c_Denlat_plus_denlon_3d

c_Denlat_plus_denlon_3d=None


def sum_over_l_for_omega(omega_i):
       
    kDensity1 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity2 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity3 = np.zeros((nr,nlat,nlon),dtype=complex)
    kDensity4 = np.zeros((nr,nlat,nlon),dtype=complex)
    
    kSpeed1 = np.zeros((nr,nlat,nlon),dtype=complex)
    kSpeed2 = np.zeros((nr,nlat,nlon),dtype=complex)
    
    for ell in xrange(0,10):
        
        filename = 'ell_{0:d}_omega_{1:d}'.format(ell, omega_i)
        matfile = os.path.join(path,filename) ; 
        g_data = loadmat(matfile) ;
    
        g_data['xisrc_denkernel']=np.squeeze(g_data['xisrc_denkernel'])
        g_data['psrc_denkernel']=np.squeeze(g_data['psrc_denkernel'])
        g_data['xircv']=np.squeeze(g_data['xircv'])
        g_data['prcv']=np.squeeze(g_data['prcv'])
        g_data['xisrc_sskernel']=np.squeeze(g_data['xisrc_sskernel'])
        g_data['xircv_sskernel']=np.squeeze(g_data['xircv_sskernel'])
        g_data['psrc_sskernel']=np.squeeze(g_data['psrc_sskernel'])
        g_data['prcv_sskernel']=np.squeeze(g_data['prcv_sskernel'])

        norm=(2*ell+1)/(4*pi)
        
        
        kstart=time.time()
        
        tensorproduct.outer_and_add_density(
                        g_data['xisrc_denkernel'],LP_src_deriv[0,ell],kDensity1,
                        g_data['xircv'],LP_rcv_deriv[0,ell],kDensity2,
                        g_data['psrc_denkernel'],LP_src_deriv[1,ell],kDensity3,
                        g_data['prcv'],LP_rcv_deriv[1,ell],kDensity4,
                        norm
                        )

        
        tensorproduct.outer_and_add_speed(
                        g_data['xisrc_sskernel'],LP_src_deriv[0,ell],
                        g_data['psrc_sskernel'],LP_src_deriv[2,ell],kSpeed1,
                        g_data['xircv_sskernel'],LP_rcv_deriv[0,ell],
                        g_data['prcv_sskernel'],LP_rcv_deriv[2,ell],kSpeed2,
                        norm
                        )


        print "kernels computed in",time.time()-kstart
    
    return kDensity1*kDensity2+kDensity3*kDensity4 , kSpeed1*kSpeed2 
    
class kernels:
    def __init__(self):
        self.kSpeed=np.zeros((nr,nlat,nlon),dtype=complex)
        self.kDensity=np.zeros((nr,nlat,nlon),dtype=complex)
        self.lock = thread.allocate_lock()
    def add(self,result):
        self.lock.acquire()
        self.kDensity+=result[0]
        self.kSpeed+=result[1]
        self.lock.release()
    


#~ omega_range = range(20*24*node_num+1,20*24*(node_num+1))

omega_range=[1]

start=time.time()
print "Starting loop at time",time.asctime(time.localtime())




pool=multiprocessing.Pool(1)

k=kernels()
for omega_i in omega_range:
    pool.apply_async(sum_over_l_for_omega,(omega_i,),callback=k.add)

pool.close()
pool.join()

end=time.time()
print "Finished computation at time",time.asctime(time.localtime())
print "Time taken",end-start,"seconds"

#~ fnameDensity = 'Density_nomega_{:d}'.format(proc_id);
#~ kDensityfile = os.path.join( directory_kD, fnameDensity);
#~ 
#~ fnameSpeed = 'SoundSpeed_nomega_{:d}'.format(proc_id);
#~ kSpeedfile = os.path.join( directory_kSS, fnameSpeed);

#~ np.savez(kDensityfile,kDensity=kDensity)
#~ np.savez(kDensityfile,kSpeed=kSpeed)
