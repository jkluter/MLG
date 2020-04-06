__all__ = ['Fit_Motion','FitFunction_Motion','FitFunction_Motion_Scandir','Jacobian_Motion',\
		'Jacobian_Motion_Scandir','Partial_Motion_no_Mass','Jac_Motion_Index' ]

from  scipy.optimize import least_squares
import numpy as np
from MLG.Microlensing.const import const_Einsteinradius
from MLG.StellarMotion import getSun,t_0,stars_position_micro
from MLG.Math import sindeg,cosdeg, unitFac
import time





def Fit_Motion(obs,num_vec, time_vec,obs_error = None, sc_scandir = None, loc_vec = None, unit_obs = 'deg', \
	unit_pos = 'deg', unit_pm = 'mas/yr', unit_px = 'mas', exact = False,**kwargs):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	obs_error = None
	if obs_error is None: 
		obs_error = np.ones(len(num_vec))
	# translate units from and to mas
	unit_pos_fac	= 	unitFac('mas',unit_pos)
	unit_pm_fac 	= 	unitFac('mas',unit_pm)
	unit_px_fac		= 	unitFac('mas',unit_px)
	unit_obs_fac 	= 	unitFac(unit_obs, 'mas')

	# calculate earth vector


	# calculate sin(ra),cos(ra),sin(dec),cos(dec)
	scsc = sindeg(obs[0,0] * unit_obs_fac / 3.6e6) , cosdeg(obs[0,0] * unit_obs_fac / 3.6e6),\
			sindeg(obs[0,1] * unit_obs_fac / 3.6e6) , cosdeg(obs[0,1] * unit_obs_fac / 3.6e6)
	if loc_vec is None:
		earth = getSun(t_0 + time_vec)
		loc_vec = np.stack([(scsc[0] * earth[0] - scsc[1] *  earth[1])/ max(scsc[3], 1e-6),\
			scsc[1] * scsc[2] *  earth[0] + scsc[0] * scsc[2] *  earth[1] - scsc[3] *  earth[2]]).T
	#switch to relative coordinates by subtracting ra0,dec0 (first observation)
	
	radec_0  = obs[0,:]  
	obs_rel = (obs - radec_0.reshape(-1,2))*unit_obs_fac

	#reorder numtimeloc_vec
	ntl_vec = [num_vec,time_vec,loc_vec]

	#setup starting value for par and par0
	#par = [0.5, ra1,dec1, 0,0,0]
	par = np.zeros(max(num_vec + 1) * 5)
	par0 = np.zeros(max(num_vec + 1) * 5)
	par0[0::5] =  radec_0[0] * unit_obs_fac
	par0[1::5] =  radec_0[1] * unit_obs_fac
	q = [np.where(num_vec == i )[0] for i in range(max(num_vec + 1))]
	index_0 = [i[0] if len(i) > 0 else -1 for i in q ]
	par[0::5] = obs_rel[index_0, 0]
	par[1::5] = obs_rel[index_0, 1]
	qq = np.array([False if len(i) > 0 else True for i in q ])

	#fitting residual function
	if sc_scandir is None:
		par_res = least_squares(FitFunction_Motion,par, xtol = 3e-16 ,jac = Jacobian_Motion, \
					args = (ntl_vec, obs_rel, obs_error, scsc[3]))
	
	else: 
		par_res = least_squares(FitFunction_Motion_Scandir,par, xtol = 3e-16 ,jac = Jacobian_Motion_Scandir, \
					args = (ntl_vec, sc_scandir, obs_rel, obs_error, scsc[3]))


	#return parameter in requested units
	par_res.x[0::5] = (par_res.x[0::5] + par0[0::5]) * unit_pos_fac
	par_res.x[1::5] = (par_res.x[1::5] + par0[1::5]) * unit_pos_fac
	par_res.x[2::5] = (par_res.x[2::5] + par0[2::5]) * unit_pm_fac
	par_res.x[3::5] = (par_res.x[3::5] + par0[3::5]) * unit_pm_fac
	par_res.x[4::5] = (par_res.x[4::5] + par0[4::5]) * unit_px_fac
	if qq.any():
		par_res.x[1+5*np.where(qq)[0]] = None
		par_res.x[2+5*np.where(qq)[0]] = None
		par_res.x[3+5*np.where(qq)[0]] = None
		par_res.x[4+5*np.where(qq)[0]] = None
		par_res.x[5+5*np.where(qq)[0]] = None
	return par_res

def FitFunction_Motion(par, numtimeloc_vec, obs, obs_error, cd):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	expect coordinates in mas
	par = parameter:  5 parameter for each star
	numtimeloc_vec = explanation for obs data: Which star at wich epoch from wich position
	obs = observational imput
	cd = cos(dec)
	'''

	c2 = const_Einsteinradius * const_Einsteinradius
	if cd == 0: cd = 1e-6 #avoid devision by 0
	cd_vec = np.array([cd,1])

	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	#order parameter
	astrometry =par.reshape((-1,5))
	radec = astrometry[:,0:2]
	pm_radec = astrometry[:,2:4]
	px = astrometry[:,4]


	# setup  parameter for each observation point
	radec_vec = radec[num_vec]
	pm_vec = pm_radec[num_vec]
	px_vec = px[num_vec].reshape((-1,1))
	radec_lens_vec = radec[[np.zeros(len(num_vec), int)]]
	pm_lens_vec = pm_radec[[np.zeros(len(num_vec), int)]]
	px_lens_vec = px[[np.zeros(len(num_vec), int)]].reshape((-1,1))

	res = np.sqrt(np.sum(np.square(\
				((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
 				- obs) * cd_vec ),axis= 1) / obs_error ** 2).reshape((-1,1))

	return res 

def FitFunction_Motion_Scandir(par, numtimeloc_vec, sc_scandir, obs, obs_error, cd):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	expect coordinates in mas
	par = parameter: 5 parameter for each star
	numtimeloc_vec = explanation for obs data: Which star at wich epoch from wich position
	obs = observational imput
	cd = cos(dec)
	'''
	c2 = const_Einsteinradius * const_Einsteinradius
	if cd == 0: cd = 1e-6 #avoid devision by 0
	cd_vec = np.array([cd,1])

	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	#order parameter
	astrometry = par.reshape((-1,5))
	radec = astrometry[:,0:2]
	pm_radec = astrometry[:,2:4]
	px = astrometry[:,4]


	# setup  parameter for each observation point
	radec_vec = radec[num_vec]
	pm_vec = pm_radec[num_vec]
	px_vec = px[num_vec].reshape((-1,1))
	radec_lens_vec = radec[[np.zeros(len(num_vec), int)]]
	pm_lens_vec = pm_radec[[np.zeros(len(num_vec), int)]]
	px_lens_vec = px[[np.zeros(len(num_vec), int)]].reshape((-1,1))

	res = np.sum(sc_scandir \
				* (((radec_vec + pm_vec / cd_vec * time_vec  + px_vec * loc_vec) \
				- obs) * cd_vec),axis= 1) / obs_error
	return res 

def Jacobian_Motion(par, numtimeloc_vec, obs, obs_error, cd):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	Determine the jacobian matrix for the  Residual FitFunction.
	input:
	par = fitparameter
	numtimeloc_vec = explanation for obs data: Which star at wich epoch from wich position:
	obs = observed/simulated data
	cd = cos(dec) for the field
	returns: Jacobian matrix   dFi/dxj
	calculation:   dFi/dxj = 1/Fi*2(delta_ra_i*d_delta_ra_i/dxj+delta_dec_i*d_delta_ra_i/dxj)
	'''
	# setup variable and constants 
	if cd == 0: cd = 1e-6
	cd_vec = np.array([cd,1])
	c2 = const_Einsteinradius * const_Einsteinradius
	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	# order parameters
	astrometry = par.reshape((-1,5))
	radec = astrometry[:,0:2]
	pm_radec = astrometry[:,2:4]
	px = astrometry[:,4]

	# setup  parameter for each observation point
	radec_vec = radec[num_vec]
	pm_vec = 	pm_radec[num_vec]
	px_vec = 	px[num_vec].reshape((-1,1))
	radec_lens_vec=radec[[np.zeros(len(num_vec), int)]]
	pm_lens_vec=pm_radec[[np.zeros(len(num_vec), int)]]
	px_lens_vec=px[[np.zeros(len(num_vec), int)]].reshape((-1,1))

	# calculat partial derivertivs d_delta_ra_i/dxj,d_delta_dec_i/dxj for Motion and microlensing 
	# cross terms dF(star_i)/ dpar(star_j) for (i != j and j != 0) are equal 0
	Motion = Partial_Motion_no_Mass(time_vec, loc_vec, cd,exact)
	
		# calculate F_i (see FitFunction)
	f_of_x = np.sqrt(np.sum(np.square(\
			((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
			- obs) * cd_vec ),axis= 1) / obs_error ** 2).reshape((-1,1))

		# calculate delta_ra_i,delta_dec_i
	a_d = (((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec)\
			- obs) * [cd,1]).reshape((-1,1,2)) / obs_error	



	# calculate dFi/dxj
	partial = (1 / np.maximum(f_of_x,1e-6) * (np.sum(a_d * (Motion) , axis = 2))).reshape(-1)

	# reshape matrix to include all cross terms
	jac = np.zeros(len(par) * len(num_vec))
	ind_par, ind_jac = Jac_Motion_Index(num_vec, len(par))
	jac[ind_jac] = partial[ind_par]

	return jac.reshape(-1, len(par))

def Jacobian_Motion_Scandir(par, numtimeloc_vec, sc_scandir, obs, obs_error, cd ):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	Determine the jacobian matrix for the  Residual FitFunction.
	input:
	par = fitparameter
	numtimeloc_vec = explanation for obs data: Which star at wich epoch from wich position:
	obs = observed/simulated data
	cd = cos(dec) for the field
	returns: Jacobian matrix   dFi/dxj
	calculation:   dFi/dxj = 1/Fi*(delta_ra_i*d_delta_ra_i/dxj+delta_dec_i*d_delta_dec_i/dxj)
	'''

	# setup variable and constants 


	if cd == 0: cd = 1e-6
	cd_vec = np.array([cd,1])
	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	# order parameters
	astrometry = par.reshape((-1,5))
	radec = astrometry[:,0:2]
	pm_radec = astrometry[:,2:4]
	px = astrometry[:,4]

	# setup  parameter for each observation point
	radec_vec = radec[num_vec]
	pm_vec = pm_radec[num_vec]
	px_vec = px[num_vec].reshape((-1,1))

	# calculat partial derivertivs d_delta_ra_i/dxj,d_delta_dec_i/dxj for Motion and microlensing 
	# cross terms dF(star_i)/ dpar(star_j) for (i != j and j == 0) are equal 0
	Motion = Partial_Motion_no_Mass(time_vec, loc_vec, cd)

	# calculate dFi/dxj
	partial = np.sum(sc_scandir.reshape(-1,1,2) * Motion / obs_error.reshape(-1,1,1), axis = 2).reshape(-1)

	# reshape matrix to include all cross terms

	jac = np.zeros(len(par) * len(num_vec))
	ind_par, ind_jac = Jac_Motion_Index(num_vec, len(par))
	jac[ind_jac] = partial[ind_par]
	if np.isnan(jac.any()): print(2)

	return jac.reshape(-1, len(par))

def Partial_Motion_no_Mass(time_vec, loc_vec, cd):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	calculat partial derivertiv for the Motion term
	ra_dec = ra_dec_0+pm * t + parallax * loc_vec  
	use t,loc_vec is constant for each observation
	first 6 zeros entries are for cross terms with the lens. 
	'''
	one = np.ones(len(time_vec[:,0]))
	zero = np.zeros(len(time_vec[:,0]))	
	alpha = np.array([zero, zero, zero, zero, zero, one * cd, zero, time_vec[:,0], zero, loc_vec[:,0] * cd]) 
	delta = np.array([zero, zero, zero, zero, zero, zero, one, zero, time_vec[:,0], loc_vec[:,1]]) 
	return np.array([alpha, delta]).T

def Jac_Motion_Index(num_vec, len_par):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	#determine which indices in the jacobian are not equal 0 

	i_par= np.arange(len(num_vec)).reshape(-1,1) * 10
	i_jac= np.arange(len(num_vec)).reshape(-1,1) * len_par
	e0 = np.where(num_vec == 0)
	ne0 = np.where(num_vec != 0)
	ind_jac0 = i_jac[e0] + np.arange(0,5) + 5 * num_vec[e0].reshape(-1,1)
	ind_par0 = i_par[e0] + np.arange(5,10)
	ind_jac1 = i_jac[ne0] + np.concatenate((np.arange(0,5) * np.ones((len(ne0[0]),1), int),\
				np.arange(0,5) + 5 * num_vec[ne0].reshape(-1,1)), axis = 1) 
	ind_par1 = i_par[ne0] + np.arange(10)
	ind_jac = np.concatenate((ind_jac0,ind_jac1), axis= None)
	ind_par = np.concatenate((ind_par0,ind_par1), axis= None)
	
	return ind_par,ind_jac

