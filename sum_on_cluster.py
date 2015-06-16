from __future__ import division
import numpy as np
import os
import time
t=time.time()

directory_kD = '/scratch/samrat/kernel/2015-06-16-kd-l-256src_lon_piby6_rcv_lon_piby3-nu-1.5mHzto4.7mHz-div-1580'
directory_kSS = '/scratch/samrat/kernel/2015-06-16-kss-l-256src_lon_piby6_rcv_lon_piby3-nu-1.5mHzto4.7mHz-div-1580'
nnode=15
ppn=16
nr=133
nlat=272
nlon=2*nlat
k_density=np.zeros((nr,nlat,1),dtype=float)
k_soundspeed=np.zeros((nr,nlat,1),dtype=float)

for coreid in xrange(ppn):
	density_kernel_sum=np.zeros((nr,nlat,nlon//ppn),dtype=float)
	for nodeid in xrange(nnode):
		t1=time.time()
		filename_Density = 'Density-nodeid-{:d}-procid-{:d}.npz'.format(nodeid,coreid)
		kDensityfile = os.path.join( directory_kD, filename_Density)		
		dkernel=np.load(kDensityfile)
		density_kernel_sum+=dkernel['density_kernel']
		dkernel=None

		print 'done for density for nodeid and coreid ',nodeid,coreid,'in time',(time.time()-t1),'sec'
		

	k_density=np.dstack((k_density,density_kernel_sum))
k_density=np.delete(k_density,0,2)
np.savez(os.path.join(directory_kD,'kernel_d'),kernel_d=k_density)
k_density=None
for coreid in xrange(ppn):
	ss_kernel_sum=np.zeros((nr,nlat,nlon//ppn),dtype=float)
	for nodeid in xrange(nnode):  
		t2=time.time()              
		filename_SoundSpeed = 'SoundSpeed-nodeid-{:d}-procid-{:d}.npz'.format(nodeid,coreid)
		kSoundSpeedfile = os.path.join( directory_kSS, filename_SoundSpeed)       
		Skernel=np.load(kSoundSpeedfile)
		ss_kernel_sum+=Skernel['sspeed_kernel']
		Skernel=None
		print 'Now loading soundspeed for nodeid and coreid ',nodeid,coreid,'in time',(time.time()-t2),'sec'
	k_soundspeed=np.dstack((k_soundspeed,ss_kernel_sum))
k_soundspeed=np.delete(k_soundspeed,0,2)
np.savez(os.path.join(directory_kSS,'kernel_s'),kernel_s=k_soundspeed)
k_soundspeed=None
print 'total time taken is ',time.time()-t,'sec'
