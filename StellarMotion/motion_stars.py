import numpy as np
from astropy.coordinates import get_sun as get_sun_astropy
from astropy.time import Time
from MLG.Microlensing.const import *
from MLG.Math import sindeg,cosdeg,sinmas,cosmas, unitFac
__all__ = ['t_0','getSun','stars_position_micro','stars_position']

t_0 = 2015.5 #GaiaDR2 epoch

def getSun(T):
	#returns the position of Lagrange Point 2
	return -1.01 * get_sun_astropy(Time(T, format = 'jyear')).cartesian.xyz.to('AU').value



def stars_position_micro(par,num_vec,time_vec,loc_vec = None,scsc = None, \
			out_unit = 'mas',unit_pos = 'deg', unit_pm = 'mas/yr' , unit_px = 'mas', \
			ml = True, pml = False, exact = True, **kwargs):
	'''----------------------------------------------------------------
	Calculate the postion of a list of stars at given times. 
	-------------------------------------------------------------------
		Input:
			par: 	  parameter of the lens and surcees [mass, ra,dec,pmra,pmdec,px, ra_source,dec_source,...]
					  or list of parameters [[mass, ra,dec,pmra,pmdec,px, ra_source,dec_source,...],...]
					  (for varying the input parameters)
			num_vec:  identifier does the observation coresponds to the lens or source
			time_vec: list of dates for the observvation
			loc_vec:  precalculate position of Gaia respect to the sun
			scsc: pre calculate sin(ra),cos(ra), sin(dec), cos(dec)
			out_unit: unit of the output
			unit_pos: unit of the position
			unit_pm:  unit of the propermotion
			unit_px:  unit of the parallax
			ml:		  Include microlensing in the calculation
			pml:	  Calculate the photometric magnification
			exact: 	  use aproximation or exact formular for pos_+  	
	-------------------------------------------------------------------
		Output:
			radec_obs:  array or matrix of [ra, dec] for each time  
			A: 			if pml: magnification for each of the observations.
	----------------------------------------------------------------'''
	
	#----------------------------------------------------------------
	#calculate scalfactors 
	unit_pos_fac	= 	unitFac(unit_pos,'mas')
	unit_pm_fac 	= 	unitFac(unit_pm, 'mas')
	unit_px_fac		= 	unitFac(unit_px, 'mas')
	unit_out_fac 	= 	unitFac('mas', out_unit)
	#----------------------------------------------------------------

	if isinstance(par[0], (int,float)): 	#check if only one set of parameters is given:
		#----------------------------------------------------------------
		# reorder parameters  
		par = np.array(par)	
		mass = par[0]

		if len(par)%5 == 1: # check if mass is given in the list of parameters 
			astrometry =par[1:].reshape((-1,5))
		else: 
			astrometry =par.reshape((-1,5))	

		radec = astrometry[:,0:2] * unit_pos_fac
		pm_radec = astrometry[:,2:4] * unit_pm_fac
		px = astrometry[:,4] * unit_px_fac	
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		#calculate scsc and loc_vec, if not given 
		if scsc is None: 
			scsc = sinmas(radec[0,0]*unit_pos_fac),cosmas(radec[0,0]),sinmas(radec[0,1]),\
				min(cosmas(radec[0,1]),1e-16),
		if loc_vec is None:
			earth = getSun(t_0 + time_vec).T
			loc_vec = np.stack([(scsc[0] * earth[:,0] - scsc[1] *  earth[:,1])\
					/ max(scsc[3],1e-6),\
					scsc[1] * scsc[2] *  earth[:,0] + scsc[0] * scsc[2]\
					*  earth[:,1] - scsc[3] *  earth[:,2]]).T.reshape(-1,1,2)
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		# reorder parameters  (cont.)
		radec_vec = radec[num_vec]	
		pm_vec = 	pm_radec[num_vec] 
		px_vec = 	px[num_vec].reshape((-1,1))
		radec_lens_vec = radec[np.zeros(len(num_vec),int)] 	# position of lens for all epochs  (cont.)
		pm_lens_vec = pm_radec[np.zeros(len(num_vec),int)]	# position of lens for all epochs  (cont.)
		px_lens_vec = px[[np.zeros(len(num_vec),int)]].reshape((-1,1)) # position of lens for all epochs  (cont.)
		T_vec =  time_vec.reshape((-1,1)) 
		cd = np.array([scsc[3],1])
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		# calculate unlensed positon:
		radec_cur  = radec_vec+ pm_vec/cd  * T_vec + px_vec * loc_vec
		radec_lens_cur = radec_lens_vec + pm_lens_vec/cd * T_vec + px_lens_vec  * loc_vec
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		# calculate shif and Magnification due to microlensing 
		A= 1 # default value for the Magnification  
		if ml: 
			#calculate Einstein radius squared
			ThetaE2 = const_Einsteinradius*const_Einsteinradius * mass*(px_lens_vec-px_vec).reshape((-1,1))
			
			if exact: 
				radec_obs = radec_cur * unit_out_fac + unit_out_fac * (radec_cur-radec_lens_cur) \
					*(np.sqrt(1/4. + ThetaE2 / np.maximum(
					np.sum(np.square((radec_cur-radec_lens_cur)* cd),axis = 1),1e-16).reshape((-1,1))) -1 /2.)
			else: 
				radec_obs = radec_cur * unit_out_fac + (radec_cur-radec_lens_cur) * ThetaE2 * unit_out_fac \
					/ np.maximum((np.sum(np.square((radec_cur-radec_lens_cur)* cd),\
					axis = 1).reshape((-1,1))+ 2 * ThetaE2),1e-16) 
			if pml:
				dist2 = np.sum(np.square((radec_cur-radec_lens_cur)* cd),axis = 1).reshape((-1,1))
				A = (dist2+2 * ThetaE2)/np.sqrt(np.maximum(dist2*dist2 + 4 * dist2 * ThetaE2,1e-16))

		else: radec_obs = radec_cur * unit_out_fac #without lensing
		#----------------------------------------------------------------
		if pml:
			return radec_obs, A
		else:
			return radec_obs
	
	else: #If a matrix of parameters is given 
		#----------------------------------------------------------------
		# reorder parameters  
		par = np.array(par)
		mass = par[:,0]

		if par.shape[1]%5 == 1: 
			astrometry =par[:,1:].reshape(len(par),-1,5)
		else: 
			astrometry =par.reshape(len(par),-1,5)

		radec = astrometry[:,:,0:2] * unit_pos_fac
		pm_radec = astrometry[:,:,2:4] * unit_pm_fac
		px = astrometry[:,:,4] * unit_px_fac
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		#calculate scsc and loc_vec, if not given 
		if scsc is None:
			scsc = np.vstack([sinmas(radec[:,0,0]*unit_pos_fac),cosmas(radec[:,0,0]), \
					sinmas(radec[:,0,1]),np.minimum(cosmas(radec[:,0,1]),1e-16)]).T
		if loc_vec is None:
			earth = getSun(t_0 + time_vec).T
			loc = np.stack([(scsc[:,0] * earth[:,:,0] - scsc[:,1] *  earth[:,:,1])\
					/ np.maximum(scsc[:,3],1e-6),\
					scsc[:,1] * scsc[:,2] *  earth[:,:,0] + scsc[:,0] * scsc[:,2]\
					*  earth[:,:,1] - scsc[:,3] *  earth[:,:,2]]).reshape(len(par),-1,1,2)
		#----------------------------------------------------------------
		
		#----------------------------------------------------------------
		# reorder parameters  (cont.)
		radec_vec = radec[:,num_vec]
		pm_vec = 	pm_radec[:,num_vec]
		px_vec = 	px[:,num_vec].reshape((len(par),-1,1))
		radec_lens_vec = radec[:,np.zeros(len(num_vec),int)]
		pm_lens_vec = pm_radec[:,np.zeros(len(num_vec),int)]
		px_lens_vec = px[:,[np.zeros(len(num_vec),int)]].reshape((len(par),-1,1))
		T_vec =  np.repeat(time_vec.reshape((1,-1,1)), len(par),axis = 0)
		cd = np.array([scsc[3],1])
		#----------------------------------------------------------------

		#----------------------------------------------------------------
		# calculate unlensed positon:
		radec_cur  = radec_vec+ pm_vec/cd  * T_vec + px_vec * loc_vec
		radec_lens_cur = radec_lens_vec + pm_lens_vec/cd * T_vec + px_lens_vec  * loc_vec
		#----------------------------------------------------------------
		
		#----------------------------------------------------------------
		# calculate shif and Magnification due to microlensing 
		A= 1 # default value for the Magnification  
		if ml: 
			#calculate Einstein radius squared
			ThetaE2 = const_Einsteinradius*const_Einsteinradius * mass.reshape(-1,1,1)\
			 			* (px_lens_vec-px_vec).reshape((len(par),-1,1))
			if exact: 
				radec_obs = radec_cur * unit_out_fac + unit_out_fac * (radec_cur-radec_lens_cur) \
					*(np.sqrt(1/4. + ThetaE2 / np.maximum(
					np.sum(np.square((radec_cur-radec_lens_cur)* cd),axis = 2),1e-16).reshape((len(par),-1,1))) -1 /2.)
			else: 

				radec_obs = radec_cur * unit_out_fac + (radec_cur-radec_lens_cur) * ThetaE2 * unit_out_fac \
					/ np.maximum((np.sum(np.square((radec_cur-radec_lens_cur)* cd),\
					axis = 2).reshape((len(par),-1,1))+ 2 * ThetaE2),1e-16)
			if pml:
				dist2 = np.sum(np.square((radec_cur-radec_lens_cur)* cd),axis = 2).reshape((len(par),-1,1))
				A = (dist2+2 * ThetaE2)/np.sqrt(np.maximum(dist2*dist2 + 4 * dist2 * ThetaE2,1e-16))
		else: radec_obs = radec_cur * unit_out_fac
		#----------------------------------------------------------------
		if pml:
			return radec_obs, A 
		else:
			return radec_obs 

def stars_position(par,num_vec,time_vec,loc_vec = None,scsc = None,
			out_unit = 'mas',unit_pos = 'deg', unit_pm = 'mas/yr' , unit_px = 'mas',**kwargs):
	#shortcut to determine the position of the star without microlensing
	return stars_position_micro(par,num_vec,time_vec, ml = False, loc_vec = None,scsc = None,
			out_unit = 'mas',unit_pos = 'deg', unit_pm = 'mas/yr' , unit_px = 'mas', **kwargs)



