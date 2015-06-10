#changes made
#constant damping 50 micro Hz.
#datas are saved on a reduced grid of 133 points
from __future__ import division
from __future__ import print_function
import numpy as np
import mi_cython as mc
import os
from scipy import interpolate
import errno
import time
from time import gmtime, strftime
nnode=15
ppn=16
ellmax =480
numax = 5.5
numin = 0.5
damping_by_resolution_factor=10.0
damping_const=50e-6 #50microHz 

pi = np.pi
exp = np.exp
sqrt = np.sqrt
codedir='/home/samrat/SoKer'#os.getcwd()
nproc=ppn*nnode
procid = int(os.environ['PBS_VNODENUM'])

sampling = 2 
rsun = 6.9598e10
diml = rsun
dimc = 1e6
dimrho = 1e-3
width = (50e5)/rsun
#center_src = 1. -(200e5)/rsun
#center_rcv = 1. -(200e5)/rsun

def derivfd(y):
      dy = np.zeros(np.size(y))
      dy = 0.5 * ( np.roll(y,-1) - np.roll(y,1))
      dy[0] = y[1]-y[0]
      dy[np.size(y)-1] = y[np.size(y)-1]-y[np.size(y)-2]      
      return dy
ndiv=int(damping_by_resolution_factor*(numax-numin)*1e-3/damping_const)
ndiv = int(ndiv-ndiv%10)
frequency_resolution=(numax-numin)*1e-3/ndiv

print ('frequency resolution-',frequency_resolution,'muHz- damping',damping_const,'muHz- divisions',ndiv)

nu_real = np.linspace( numin, numax, ndiv) *1e-3
#~ dampingfile=os.path.join(codedir,"m585q.4816")
#~ dmpnu = np.loadtxt(dampingfile,usecols= [2])*1e-6
#~ fwhm = np.loadtxt(dampingfile,usecols= [4])*1e-6 *1e1 

# Sort the damping array before interpolating onto our grid
#~ points = zip(dmpnu, fwhm)
#~ points = sorted(points, key=lambda point: point[0])
#~ dmpnu, fwhm = zip(*points)
#~ 
#~ damping = interpolate.interp1d(dmpnu, fwhm) 
nu = (nu_real + 1j*damping_const)*(diml/dimc) #mHz 

nu0 = 3e-3 *diml/dimc #mHz
sigma = 1e-3 *diml/dimc
dnu = derivfd(nu_real)
fnu = (-1)*(nu_real**2)*exp(-(nu_real-nu0)**2 / (2*sigma**2))

solar_model=os.path.join(codedir,"model_S(GONG)")
r_true, rho,c, gravity= np.loadtxt(solar_model, unpack=True)
r = r_true[::sampling]
nr = np.size(r)
c = c[::sampling]
c2 = c**2 / dimc**2
gravity=gravity[::sampling]
rho=rho[::sampling]/dimrho

r_reduced = r_true[::12]
center_src = r_reduced[121]
center_rcv = r_reduced[121]
print ('center_src and rcv ', center_src, center_rcv)

scaling = derivfd(r)
goverc2 = gravity/c2
hrho = derivfd(np.log(rho))/scaling
N2 = (-1)* gravity*( goverc2 + hrho)
goverc2 = goverc2 + hrho*0.5

rho0 = interpolate.interp1d(r,rho)(r_reduced)   #rho now approximated on the reduced grid of radius
c2_reduced = interpolate.interp1d(r,c2)(r_reduced) #c2 interpolated on reduced grid of radius
scaling_reduced = interpolate.interp1d(r,scaling)(r_reduced)

rho = 1.0

src = exp((-1)*(r-center_src)**2 / (2.0*width**2))
rcv = exp((-1)*(r-center_rcv)**2 / (2.0*width**2))
src_2x = np.append(np.zeros(nr-1), src)
src_2x = np.delete(src_2x, (2*nr-2))
rcv_2x = np.append(np.zeros(nr-1), rcv)
rcv_2x = np.delete(rcv_2x, (2*nr-2))

save_path = '/scratch/samrat/kernel/'+strftime("%Y-%m-%d", gmtime())+'|ell-'+str(ellmax)+'-frequency-'+str(numin)+'mHzto'+str(numax)+'mHz-divisions-'+str(ndiv)
if not os.path.isdir(save_path):
	try:
		os.makedirs(save_path)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise e
		pass
if os.path.isfile('temp')==False:
	info=open('temp','w')
	print(save_path,'\t',numin,'\t',numax,'\t',ndiv,'\t',ellmax,file=info)
	info.close()
os.chdir(save_path)
if os.path.isfile('readme')==False:
	info=open('readme','w')
	print('numin:',numin,'\n','numax:',numax,'\n','frequency divisions:',ndiv,'\n','ellmax:',ellmax,'\n',
			'r grid computed:',len(r),'\n','r grid saved:',len(r_reduced),file=info)
	info.close()
omegalist=np.arange(procid,ndiv,nproc)  #omega list modified to run the code faster
print ('Procid-',procid,' frequencies-',omegalist)

for omegai in omegalist:
	t=time.time()
	omega = nu [ omegai] * 2*pi
	rhoomega2 = rho*(omega**2 - N2)
	xisrc, psrc, xisrcatrcv, xisrcatsrc, psrcatrcv, xircv, prcv, xisrc_denkernel ,xisrc_sskernel, psrc_denkernel, psrc_sskernel,xircv_sskernel, prcv_sskernel =  ([] for i in range(13))
	filename = "omega-"+str(omegai).zfill(4)
    
	for ell in xrange(ellmax+1):
        #print ell
		
		sl2 = ell*(ell+1)*c2/r**2
		oneoverrhoc2 = (sl2 /omega**2 -1)/ (rho*c2)

		my_matrix = mc.genmatrix( goverc2, oneoverrhoc2, rhoomega2, r)
		my_matrix = np.delete(my_matrix, [0,2*nr-1],0)
		my_matrix = np.delete(my_matrix, [0,2*nr-1],1)

		solsrc = np.linalg.solve(my_matrix, src_2x)
		xis , ps = np.split(solsrc,2)
		xis = np.append(0,xis)
		ps = np.append(ps, 0)
		#-------------------------------------------------------------------
		xis=interpolate.interp1d(r,xis)(r_reduced)      #
		ps=interpolate.interp1d(r,ps)(r_reduced)        #
		#-------------------------------------------------------------------        
		xisrc.append( xis/sqrt(rho0))
		psrc.append( ps*sqrt(rho0))
		solrcv = np.linalg.solve(my_matrix, rcv_2x)
		xir , pr = np.split(solrcv,2)
		xir = np.append(0,xir)
		pr = np.append(pr, 0)
		#--------------------------------------------------------------------
		xir=interpolate.interp1d(r,xir)(r_reduced)      #
		pr=interpolate.interp1d(r,pr)(r_reduced)        #
		#--------------------------------------------------------------------
		xircv.append( xir/sqrt(rho0))
		prcv.append( pr*sqrt(rho0))
		#---------------------------------------------------------------------
		xisrcatrcv.append(interpolate.interp1d(r_reduced, xisrc[ell])(center_rcv))      #
		xisrcatsrc.append(interpolate.interp1d(r_reduced, xisrc[ell])(center_src))      #
		psrcatrcv.append(interpolate.interp1d(r_reduced, psrc[ell])(center_rcv))        #
		#----------------------------------------------------------------------
		xisrckernel = ((2*pi)**3 *fnu[omegai]*dnu[omegai]) * xisrc[ell]
		xisrc_denkernel.append( xisrckernel * rho0)

				
		dxisrc_dr = derivfd( r_reduced**2. * xisrckernel )/( r_reduced**2 * scaling_reduced )           #
		xisrc_sskernel.append( c2_reduced*dxisrc_dr )                                           #

		psrckernel = ((2*pi)**3 *fnu[omegai]* dnu[omegai]) * psrc[ell]
		psrc_denkernel.append((psrckernel / (r_reduced **2 * rho0 ))/ omega**4)                 #
		psrc_sskernel.append((psrckernel * c2_reduced / (r_reduced**2 *rho0)) / omega**2)       #
		#---------------------------------------------------------------------

		xircv_sskernel.append( derivfd( r_reduced**2 *xircv[ell] ) / (r_reduced**2 * scaling_reduced) )  #
		prcvkernel = ((2*pi)**3 *fnu[omegai]* dnu[omegai] )* prcv[ell]                  #
		prcv_sskernel.append((prcv[ell] / (r_reduced**2. *rho0) )/ omega**2)           #
	np.savez(filename, r=r_reduced, xisrc=xisrc,xisrcatrcv =xisrcatrcv, xisrcatsrc=xisrcatsrc,
					psrc=psrc, psrcatrcv=psrcatrcv, xircv =xircv, prcv=prcv, xisrc_denkernel=xisrc_denkernel,
					xisrc_sskernel=xisrc_sskernel, psrc_denkernel=psrc_denkernel, xircv_sskernel=xircv_sskernel,
					psrc_sskernel=psrc_sskernel, prcv_sskernel = prcv_sskernel)
	procinfo=open('procid-'+str(procid),'w')
	print('omega-list:',omegalist,'\n','finished for omega:',omegai,'\n','time:',(time.time()-t),'\n','running computation now for omega:',(omegai+1),'\n',file=procinfo)
	procinfo.close()
