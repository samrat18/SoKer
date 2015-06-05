import numpy as np

directory_kSS = '/scratch/samrat/kernel/greens/sound_speed_individual/'
directory_kD = '/scratch/samrat/kernel/greens/density_individual/'


nnode=10
ppn=24
nr=798
nlat=480
nlon=960

k_density=[]
k_soundspeed=[]
for coreid in xrange(ppn):
    
    density_kernel_sum=np.zeros((nr,nlat,nlon//ppn),dtype=float)
    ss_kernel_sum=np.zeros((nr,nlat,nlon//ppn),dtype=float)
    
    for nodeid in xrange(nnode):
        
        filename_Density = 'Density-nodeid-{:d}-procid-{:d}'.format(coreid,nodeid)
        kDensityfile = os.path.join( directory_kD, filename_Density)
        filename_SoundSpeed = 'SoundSpeed-nodeid-{:d}-procid-{:d}'.format(coreid,nodeid)
        kSoundSpeedfile = os.path.join( directory_kSS, filename_SoundSpeed)
        
        print 'Now loading for nodeid and coreid ',nodeid,coreid

        dkernel=np.load(filename_Density)
        density_kernel_sum+=dkernel['density_kernel']
        dkernel=None
        Skernel=np.load(filename_SoundSpeed)
        ss_kernel_sum+=['sspeed_kernel']
        Skernel=None
    k_density.append(density_kernel_sum)
    k_soundspeed.append(ss_kernel_sum)
