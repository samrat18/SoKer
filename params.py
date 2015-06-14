from __future__ import division
import numpy as np
import os

def params():
	codedir=os.getcwd()+'/SoKer'
	nnode=15
	ppn=16
	ellmax =480
	numin = 1.5
	numax = 4.7
	damping_by_resolution_factor=4.0
	reduced_sampling=12
	
	
	dampingfile=os.path.join(codedir,"m585q.4816")
	dmpnu = np.loadtxt(dampingfile,usecols= [2])*1e-6
	fwhm = np.loadtxt(dampingfile,usecols= [4])*1e-6 *50e1
	
	sampling = 2 
	rsun = 6.9598e10
	width = (50e5)/rsun
	#center_src = 1. -(200e5)/rsun
	#center_rcv = 1. -(200e5)/rsun
	return codedir, nnode, ppn, ellmax, numin, numax, damping_by_resolution_factor,reduced_sampling, dmpnu,fwhm,sampling,width
	
if __name__=='__params__':
    sys.exit(params())
