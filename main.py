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
import params

pi = np.pi
exp = np.exp
sqrt = np.sqrt

codedir, nnode, ppn, ellmax, numin, numax, dmp_by_res_factor,reduced_sampling, dmpnu,fwhm,sampling,width=params.params()

procid = int(os.environ['PBS_VNODENUM'])
nproc=ppn*nnode

damping_const=abs(fwhm).min()
ndiv=(dmp_by_res_factor*(numax-numin)*1e-3/damping_const)
ndiv=int(ndiv-ndiv%10)
frequency_resolution=(numax-numin)*1e-3/ndiv
print ('frequency resolution-',frequency_resolution,'muHz- damping',damping_const,'muHz',' -ndiv-',ndiv)

points = zip(dmpnu, fwhm)
points = sorted(points, key=lambda point: point[0])
dmpnu, fwhm = zip(*points)
damping = interpolate.interp1d(dmpnu, fwhm)

def derivfd(y):
      dy = np.zeros(np.size(y))
      dy = 0.5 * ( np.roll(y,-1) - np.roll(y,1))
      dy[0] = y[1]-y[0]
      dy[np.size(y)-1] = y[np.size(y)-1]-y[np.size(y)-2]      
      return dy

rsun = 6.9598e10
diml = rsun
dimc = 1e6
dimrho = 1e-3

nu_real = np.linspace( numin, numax, ndiv) *1e-3
nu = (nu_real + 1j*damping_const)*(diml/dimc) 

nu0 = 3e-3 *diml/dimc
sigma = 1e-3 *diml/dimc
dnu = derivfd(nu_real)
fnu = (-1)*exp( -(nu_real-nu0)**2/(2*sigma**2) )
fdnu=dnu*fnu*(diml/dimc)

solar_model=os.path.join(codedir,"model_S(GONG)")
r_true, rho,c, gravity= np.loadtxt(solar_model, unpack=True)
r = r_true[::sampling]
nr = np.size(r)
c = c[::sampling]
c2 = c**2 / dimc**2
gravity=gravity[::sampling]
rho=rho[::sampling]/dimrho

r_reduced = r_true[::reduced_sampling]
center_src = r_reduced[121]
center_rcv = r_reduced[121]

print ('center_src and rcv ', center_src, center_rcv)

scaling = derivfd(r)
goverc2 = gravity/c2
hrho = derivfd(np.log(rho))/scaling
N2 = (-1)* gravity*( goverc2 + hrho)
goverc2 = goverc2 + hrho*0.5

rho0 = interpolate.interp1d(r,rho)(r_reduced)   
c2_reduced = interpolate.interp1d(r,c2)(r_reduced)
gravity_reduced = interpolate.interp1d(r,gravity)(r_reduced) 
scaling_reduced = interpolate.interp1d(r,scaling)(r_reduced)

rho = 1.0

src = exp((-1)*(r-center_src)**2 / (2.0*width**2))
rcv = exp((-1)*(r-center_rcv)**2 / (2.0*width**2))
src_2x = np.append(np.zeros(nr-1), src)
src_2x = np.delete(src_2x, (2*nr-2))
rcv_2x = np.append(np.zeros(nr-1), rcv)
rcv_2x = np.delete(rcv_2x, (2*nr-2))

save_path = '/scratch/samrat/kernel/'+strftime("%Y-%m-%d", gmtime())+'-l-'+str(ellmax)+'-f-'+str(numin)+'to'+str(numax)+'-nnu-'+str(ndiv)
if not os.path.isdir(save_path):
	try:
		os.makedirs(save_path)
	except OSError, e:
		if e.errno != errno.EEXIST:
			raise e
		pass
if os.path.isfile(os.path.join(codedir,'temp'))==False:
	info=open('temp','w')
	print(save_path,'\t',numin,'\t',numax,'\t',ndiv,'\t',ellmax,file=info)
	info.close()
os.chdir(save_path)
if os.path.isfile('readme')==False:
	info=open('readme','w')
	print('numin:',numin,'\n','numax:',numax,'\n','frequency divisions:',ndiv,'\n','ellmax:',ellmax,'\n',
			'r grid computed:',len(r),'\n','r grid saved:',len(r_reduced),file=info)
	info.close()
omegalist=np.arange(procid,ndiv,nproc)
print ('Procid-',procid,' frequencies-',omegalist)

for omegai in omegalist:
	t=time.time()
	omega = nu [omegai] * 2*pi
	rhoomega2 = rho*(omega**2 - N2)
	xisrc, psrc, xisrcatrcv, psrcatrcv, xircv, prcv=([] for i in range(6)) 
	dknl_i_Grr,dknl_i_Ghr,dknl_i_Chirr,dknl_i_Chihr=([] for i in range(4))
	dknl_j_Grr,dknl_j_Ghr,dknl_j_Chirr,dknl_j_Chihr=([] for i in range(4))
	dknl_k_Gphir,dknl_k_Chirr,dknl_k_Chihr=([] for i in range(3))
	dknl_l_Grr,dknl_l_Ghr,dknl_l_Chirr,dknl_l_Chihr=([] for i in range(4))
	ssknl_Grr,ssknl_Ghr,ssknl_Chirr,ssknl_Chihr=([] for i in range(4))
	filename = "omega-"+str(omegai).zfill(4)
    
	for ell in xrange(ellmax+1):
		t=time.time()
		sl2 = ell*(ell+1)*c2/r**2
		oneoverrhoc2 = (sl2 /omega**2 -1)/ (rho*c2)

		my_matrix = mc.genmatrix( goverc2, oneoverrhoc2, rhoomega2, r)
		my_matrix = np.delete(my_matrix, [0,2*nr-1],0)
		my_matrix = np.delete(my_matrix, [0,2*nr-1],1)

		solsrc = np.linalg.solve(my_matrix, src_2x)
		solrcv = np.linalg.solve(my_matrix, rcv_2x)
		
		xis , ps = np.split(solsrc,2)
		xir , pr = np.split(solrcv,2)
		
		xis = np.append(0,xis)
		ps = np.append(ps, 0)
		xir = np.append(0,xir)
		pr = np.append(pr, 0)
		
		xis=(interpolate.interp1d(r,xis.real)(r_reduced))/sqrt(rho0)      
		ps=(interpolate.interp1d(r,ps.real)(r_reduced))*sqrt(rho0) 
		xir=(interpolate.interp1d(r,xir.real)(r_reduced))/sqrt(rho0)
		pr=(interpolate.interp1d(r,pr.real)(r_reduced))*sqrt(rho0)         
		
		xisrc.append(xis)
		psrc.append(ps)
		xircv.append(xir)
		prcv.append(pr)		

		#~ xisrcatrcv.append(interpolate.interp1d(r_reduced, xisrc[ell])(center_rcv)) 
		#~ psrcatrcv.append(interpolate.interp1d(r_reduced, psrc[ell])(center_rcv))       
		#~ xisrcatsrc.append(interpolate.interp1d(r_reduced, xisrc[ell])(center_src)) ################################
		xisrcatrcv.append(interpolate.interp1d(r_reduced,xis)(center_rcv)) 
		psrcatrcv.append(interpolate.interp1d(r_reduced, ps)(center_rcv))
		
		#G-> Greens function by source
		#Chi->Greens function by receiver
		rrhoOmega2=(r_reduced*rho0)*omega**2.
		Grr=xis *(fdnu[omegai]*2*pi)
		Ghr=ps/rrhoOmega2 *(fdnu[omegai]*2*pi)
		Chirr=xir
		Chihr=pr/rrhoOmega2
		
		div_Grr=derivfd(r_reduced**2.*Grr )/(r_reduced**2*scaling_reduced )
		div_Ghr=Ghr/r_reduced
		div_Chirr=derivfd(r_reduced**2.*Chirr)/(r_reduced**2*scaling_reduced )  
		div_Chihr=Chihr/r_reduced
		
		dknl_i_Grr.append(np.real(omega**2.) *Grr)
		dknl_i_Ghr.append(np.real(omega**2.) *Ghr)
		dknl_i_Chirr.append(Chirr)
		dknl_i_Chihr.append(Chihr)
		
		dknl_j_Grr.append(div_Grr*c2_reduced)
		dknl_j_Ghr.append(div_Ghr*c2_reduced)
		dknl_j_Chirr.append(div_Chirr)
		dknl_j_Chihr.append(div_Chihr)
		
		dknl_k_Gphir.append(gravity_reduced*Ghr)
		dknl_k_Chirr.append(div_Chirr)
		dknl_k_Chihr.append(div_Chihr)
		
		dknl_l_Grr.append(Grr)		
		dknl_l_Ghr.append(Ghr)
		dknl_l_Chirr.append(derivfd(gravity_reduced*Chirr)/scaling_reduced)
		dknl_l_Chihr.append(gravity_reduced*Chihr/r_reduced)
		
		ssknl_Grr.append(div_Grr*rho0)
		ssknl_Ghr.append(div_Ghr*rho0)
		ssknl_Chirr.append(div_Chirr)
		ssknl_Chihr.append(div_Chihr)
		
		print(time.time()-t)
		
		
		#~ xisrckernel = ((2*pi)**3 *fnu[omegai]*dnu[omegai]) * xis
		#~ xisrc_denkernel.append( xisrckernel * rho0)
#~ 
		#~ dxisrc_dr = derivfd( r_reduced**2. * xisrckernel )/( r_reduced**2 * scaling_reduced )           
		#~ xisrc_sskernel.append( rho0*dxisrc_dr )                                           
#~ 
		#~ psrckernel = ((2*pi)**3 *fnu[omegai]* dnu[omegai]) * ps
		#~ psrc_denkernel.append((psrckernel / (r_reduced **2 * rho0 ))/ omega**4)                 
		#~ psrc_sskernel.append((psrckernel/ r_reduced**2) / omega**2)       
		
		#~ xircv_sskernel.append( derivfd( r_reduced**2 *xir ) / (r_reduced**2 * scaling_reduced) )  
		#~ prcvkernel = ((2*pi)**3 *fnu[omegai]* dnu[omegai] )* pr                  
		#~ prcv_sskernel.append((pr / (r_reduced**2. *rho0) )/ omega**2)           
	np.savez(filename, r=r_reduced, xisrc=xisrc, psrc=psrc,xircv=xircv,prcv=prcv,
					xisrcatrcv =xisrcatrcv, psrcatrcv=psrcatrcv, 
					dknl_i_Grr=dknl_i_Grr,dknl_i_Ghr=dknl_i_Ghr,dknl_i_Chirr=dknl_i_Chirr,dknl_i_Chihr=dknl_i_Chihr,
					dknl_j_Grr=dknl_j_Grr,dknl_j_Ghr=dknl_j_Ghr,dknl_j_Chirr=dknl_j_Chirr,dknl_j_Chihr=dknl_j_Chihr,
					dknl_k_Gphir=dknl_k_Gphir,dknl_k_Chirr=dknl_k_Chirr,dknl_k_Chihr=dknl_k_Chihr,
					dknl_l_Grr=dknl_l_Grr,dknl_l_Ghr=dknl_l_Ghr,dknl_l_Chirr=dknl_l_Chirr,dknl_l_Chihr=dknl_l_Chihr,
					ssknl_Grr=ssknl_Grr,ssknl_Ghr=ssknl_Ghr,ssknl_Chirr=ssknl_Chirr,ssknl_Chihr=ssknl_Chihr)
	
	procinfo=open('procid-'+str(procid),'w')
	print('omega-list:',omegalist,'\n','finished for omega:',omegai,'\n','time:',(time.time()-t),'\n','running computation now for omega:',(omegai+1),'\n',file=procinfo)
	procinfo.close()
