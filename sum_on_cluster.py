import numpy as np
import os
import time
t=time.time()
directory_kSS = '/scratch/samrat/kernel/greens/sound_speed_individual_256_ells/'
directory_kD = '/scratch/samrat/kernel/greens/density_individual_256_ells/'
nnode=10
ppn=24
nr=798
nlat=288
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
		print 'Now loading density for nodeid and coreid ',nodeid,coreid,'in time',(time.time()-t1),'sec'
	k_density=np.dstack((k_density,density_kernel_sum))
k_density=np.delete(k_density,0,2)
np.savez(os.path.join(directory_kD,'density_K'),k_density=k_density)
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
np.savez(os.path.join(directory_kSS,'Sspeed_K'),k_soundspeed=k_soundspeed)
k_soundspeed=None
print 'total time taken is ',time.time()-t,'sec'
