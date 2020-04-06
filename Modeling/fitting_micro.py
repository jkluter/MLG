__all__ = ['plotmicro','Fit_Micro','FitFunction_Micro','FitFunction_Micro_Scandir','Jacobian_Micro',\
		'Jacobian_Micro_Scandir','Partial_Micro_no_Mass','Jac_Micro_Index' ]

from  scipy.optimize import least_squares
import numpy as np
try: import MLG.Modeling.fitting_motion as fm
except: pass 
from MLG.Microlensing.const import const_Einsteinradius
from MLG.StellarMotion import getSun,t_0, stars_position_micro
from MLG.Math import sindeg,cosdeg, unitFac
import matplotlib.pyplot as plt
import time
from MLG import imagepath, paperpath


def plotmicro(obs,par, scsc,trange=[-5,10,1000],**kwargs):
	'''------------------------------------------------------------
		Description: plot the position, of lens and source star, as well as the fitted curve 
	---------------------------------------------------------------
		Input: 
			obs: set of observations [[ra,dec],...]
			par: fitted parameter
	------------------------------------------------------------'''
	time_vec = np.linspace(trange[0],trange[1], num=trange[2], endpoint=True)
	earth = getSun(t_0 + time_vec)
	loc_vec = np.stack([(scsc[0] * earth[0] - scsc[1] *  earth[1])/ max(scsc[3],1e-6),\
			scsc[1] * scsc[2] *  earth[0] + scsc[0] * scsc[2] *  earth[1] - scsc[3] *  earth[2]]).T
	n_stars = int((len(par)-1)/5)
	plt.plot(obs[:,0],obs[:,1],'x')
	styl = kwargs.pop('styl', '-')
	linewidth = kwargs.pop('linewidth', 2)
	for i in range(n_stars):
		xx = stars_position_micro(par,np.ones(1000, int)*i,time_vec,loc_vec,scsc,**kwargs)
		yy = stars_position_micro(par,np.ones(1000, int)*i,time_vec,loc_vec,scsc,ml = True, **kwargs)
		plt.plot(xx[:,0],xx[:,1],styl,linewidth=linewidth)
		plt.plot(yy[:,0],yy[:,1],':',linewidth=linewidth)


def Fit_Micro(obs, num_vec, time_vec, obs_error = None, sc_scandir = None, loc_vec = None, unit_obs = 'deg', \
	unit_pos = 'deg', unit_pm = 'mas/yr', unit_px = 'mas', plot = False, exact = False, bounds = False, prefit_motion = False,**kwargs):
	'''------------------------------------------------------------
		Description: Fited the stelar motion to de data 
	---------------------------------------------------------------
		Input:
			obs: list of observations
			numvec: identifier of the coresponding star
			time_vec: epoch of the observations
			obs_error: measurement errors for weigting. If none use uniform weights)
			sc_scandir: sin,cos of the position angle of the scan direction. for each observation
				If None 2D residuals are used.
			loc_vec: Location of Gaia. respect to the sun. If None, calculate the position from the epoch
			unit_**: Units of the observation, position, propermotion and parallax
			exact: If False, use the center of light as aproximation
			bounds: Use bounds for the fitting process. 
			prefit_motion: Determine a priors from a fit without microlensing.
	---------------------------------------------------------------
		Output:
			par_res:	leastsquare fit results
						par_res.x = [mass, ra0_Lens, dec0_Lens, pmra_Lens,pmdec_Lens,px_Lens, ra0_Source1,.... ]
	------------------------------------------------------------'''
	if obs_error is None: 
		obs_error = np.ones(len(obs)) 
	# translate units from and to mas
	unit_pos_fac	= 	unitFac('mas',unit_pos)
	unit_pm_fac 	= 	unitFac('mas',unit_pm)
	unit_px_fac		= 	unitFac('mas',unit_px)
	unit_obs_fac 	= 	unitFac(unit_obs, 'mas')



	# calculate sin(ra),cos(ra),sin(dec),cos(dec)
	scsc = sindeg(obs[0,0] * unit_obs_fac / 3.6e6) , cosdeg(obs[0,0] * unit_obs_fac / 3.6e6),\
			sindeg(obs[0,1] * unit_obs_fac / 3.6e6) , cosdeg(obs[0,1] * unit_obs_fac / 3.6e6)
	# calculate earth vector
	if loc_vec is None:
		earth = getSun(t_0 + time_vec)
		loc_vec = np.stack([(scsc[0] * earth[0] - scsc[1] *  earth[1])/ max(scsc[3], 1e-20),\
			scsc[1] * scsc[2] *  earth[0] + scsc[0] * scsc[2] *  earth[1] - scsc[3] *  earth[2]]).T
	#switch to relative coordinates by subtracting ra0,dec0 (first observation)
	
	radec_0  = obs[0,:]  
	obs_rel = (obs - radec_0.reshape(-1,2))*unit_obs_fac

	#reorder numtimeloc_vec
	ntl_vec = [num_vec,time_vec,loc_vec]

	#setup starting value for par and par0
	#par = [0.5, ra1,dec1, 0,0,0]
	if len(np.where(num_vec == 0)[0]) == 0:
		#test source only
		par = np.zeros(max(num_vec + 1) * 5 + 1)
		par0 = np.zeros(max(num_vec + 1) * 5 + 1)
		par0[1::5] =  radec_0[0] * unit_obs_fac
		par0[2::5] =  radec_0[1] * unit_obs_fac
		par[0] = 0.5
		q = [np.where(num_vec == i )[0] for i in range(max(num_vec + 1))]
		index_0 = [i[0] if len(i) > 0 else -1 for i in q ]
		par[1::5] = obs_rel[index_0, 0]
		par[2::5] = obs_rel[index_0, 1]
		par[2] = par[7]
		par[1] = par[6]
		#calculat preset parameters 
		if sc_scandir is None:
			par_1 = least_squares(fm.FitFunction_Motion,par[1:], xtol = 3e-8 ,jac = fm.Jacobian_Motion, \
					args = (ntl_vec, obs_rel, obs_error, scsc[3]))
		else:
			par_1 = least_squares(fm.FitFunction_Motion_Scandir,par[1:], xtol = 3e-8 ,jac = fm.Jacobian_Motion_Scandir, \
					args = (ntl_vec, sc_scandir, obs_rel, obs_error, scsc[3]))
		nn = np.where(np.abs(par_1.fun) == max(np.abs(par_1.fun)))[0]
		mm = par_1.fun[nn] * obs_error[nn]
		ThetaE2 = const_Einsteinradius*const_Einsteinradius * 0.5 * 100
		deltaphi = -(ThetaE2/mm + np.sqrt(max(ThetaE2**2/mm**2 - 8*ThetaE2,0)))/2
		par[1:3] = obs_rel[nn] + deltaphi * sc_scandir[nn]
		par[5] = 100 
		par[6:] = par_1.x[5:]
		qq = np.array([False,False])
	else:	
		par = np.zeros(max(num_vec + 1) * 5 + 1)
		par0 = np.zeros(max(num_vec + 1) * 5 + 1)
		par0[1::5] =  radec_0[0] * unit_obs_fac
		par0[2::5] =  radec_0[1] * unit_obs_fac
		par[0] = 0.5	
		if prefit_motion:
			if sc_scandir is None:
				par_res = least_squares(fm.FitFunction_Motion,par[1:], xtol = 3e-8 ,jac = fm.Jacobian_Motion, \
					args = (ntl_vec, obs_rel, obs_error, scsc[3]))
			else: 
				par_res = least_squares(fm.FitFunction_Motion_Scandir,par[1:], xtol = 3e-8 ,jac = fm.Jacobian_Motion_Scandir, \
					args = (ntl_vec, sc_scandir, obs_rel, obs_error, scsc[3]))
			par0[1:] = par_res.x

		q = [np.where(num_vec == i )[0] for i in range(max(num_vec + 1))]
		index_0 = [i[0] if len(i) > 0 else -1 for i in q ]
		par[1::5] = obs_rel[index_0, 0]
		par[2::5] = obs_rel[index_0, 1]
		qq = np.array([False if len(i) > 0 else True for i in q ])

	if bounds: bounds=( -np.inf, [1000, *(np.inf*np.ones(len(par)-1))])
	else:  bounds = (-np.inf, np.inf)
	#fitting residual function
	if sc_scandir is None:
		par_res = least_squares(FitFunction_Micro,par, xtol = 3e-16 ,jac = Jacobian_Micro, bounds = bounds, \
					args = (ntl_vec, obs_rel, obs_error, scsc[3], exact))
	
	else: 
		par_res = least_squares(FitFunction_Micro_Scandir,par, xtol = 3e-16 ,jac = Jacobian_Micro_Scandir, \
					bounds = bounds, args = (ntl_vec, sc_scandir, obs_rel, obs_error, scsc[3], exact))
		r = np.random.uniform()
		if False:
			cost = lambda x, i: np.sum(FitFunction_Micro_Scandir([*par_res.x[0:i],x,*par_res.x[i+1:]],ntl_vec, sc_scandir, obs_rel, obs_error, scsc[3], exact)**2)


			xx = np.linspace(-1,1,1000)
			sfactor = [1,1,1,0.1,0.1,1,1,1,0.1,0.1,1]

			y = [np.array([cost(x*sfactor[i]+par_res.x[i], i) for x in xx]) for i in range(11)]
			par_name = ['mass, M_sun', 'ra1, mas', 'dec1, mas','pmra1, 0.1mas/yr', 'pmdec1, 0.1mas/yr','px1, mas', 'ra2', 'dec2','pmra2', 'pmdec2','px2']
			fig = plt.figure(figsize = (11,10))
			for i in range(11):
				plt.plot(xx,y[i],label = par_name[i])
			plt.xlabel('delta_par') 
			plt.ylim([0,5 * par_res.cost])
			plt.legend(fontsize = 20)
			fig.savefig(imagepath + 'Cost_function_' + str(int(r*10000))+'.png')
			plt.close(fig)

			fig = plt.figure(figsize = (11,10))
			xx = np.linspace(-1000,1000,10000)
			xx2 = np.linspace(-999.999,1000.001,10000)
			y = [np.array([cost(x*sfactor[i]+par_res.x[i], i) for x in xx]) for i in range(1)]
			y2 = [np.array([cost(x*sfactor[i]+par_res.x[i], i) for x in xx2]) for i in range(1)]
			for i in range(1):
				plt.plot(xx,y[i]-y2[i],label = par_name[i])
			fig.savefig(imagepath + 'Cost_function_' + str(int(r*10000))+'_2.png')
			plt.close(fig)



	if plot:
		plotMicro(obs_rel, par_res.x, scsc, styl = '-', linewidth=2, out_unit = 'mas', \
					unit_pos = 'mas', unit_pm = 'mas/yr', unit_px = 'mas')
	#return parameter in requested units

	par_res.x[1::5] = (par_res.x[1::5] + par0[1::5]) * unit_pos_fac
	par_res.x[2::5] = (par_res.x[2::5] + par0[2::5]) * unit_pos_fac
	par_res.x[3::5] = (par_res.x[3::5] + par0[3::5]) * unit_pm_fac
	par_res.x[4::5] = (par_res.x[4::5] + par0[4::5]) * unit_pm_fac
	par_res.x[5::5] = (par_res.x[5::5] + par0[5::5]) * unit_px_fac
	if qq.any():
		par_res.x[1+5*np.where(qq)[0]] = None
		par_res.x[2+5*np.where(qq)[0]] = None
		par_res.x[3+5*np.where(qq)[0]] = None
		par_res.x[4+5*np.where(qq)[0]] = None
		par_res.x[5+5*np.where(qq)[0]] = None
	return par_res

def FitFunction_Micro(par, numtimeloc_vec, obs, obs_error, cd, exact):
	'''------------------------------------------------------------
		Description: 
			2D Residual between observation and motion including Microlensing 
			at given epochs  

	---------------------------------------------------------------
		Input:
			par:	parameter: mass plus 5 parameter for each star
					expect coordinates in mas, and mass in Solar_masses
			numtimeloc_vec: [Identifier, Epoch, GaiaLocation] for each observation
			obs:	observational data
			cd: 	cos(delta)
			exact:	if True use exact formular, else use the center of light as aproximation

	---------------------------------------------------------------
		res: 2D Residuals between Model and observation
	------------------------------------------------------------'''

	#------------------------------------------------------------
	c2 = const_Einsteinradius * const_Einsteinradius
	if cd == 0: cd = 1e-20 #avoid devision by 0
	cd_vec = np.array([cd,1])
	#------------------------------------------------------------

	#------------------------------------------------------------
	# split numtimeloc_vec 
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]
	#------------------------------------------------------------

	#------------------------------------------------------------
	#order parameter
	par = np.array(par)
	mass = par[0]
	astrometry =par[1:].reshape((-1,5))
	radec = astrometry[:,0:2]
	pm_radec = astrometry[:,2:4]
	px = astrometry[:,4]
	#------------------------------------------------------------

	#------------------------------------------------------------
	# setup parameter for each observation point
	radec_vec = radec[num_vec]
	pm_vec = pm_radec[num_vec]
	px_vec = px[num_vec].reshape((-1,1))
	radec_lens_vec = radec[[np.zeros(len(num_vec), int)]]
	pm_lens_vec = pm_radec[[np.zeros(len(num_vec), int)]]
	px_lens_vec = px[[np.zeros(len(num_vec), int)]].reshape((-1,1))
	#------------------------------------------------------------

	if exact:
		#------------------------------------------------------------
		#calculate residual:
		res = np.sqrt(np.sum(np.square(\
				((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
				+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
				+ (px_vec - px_lens_vec) * loc_vec) \
				* (np.sqrt(np.maximum(1/4. + c2 * mass * (px_lens_vec - px_vec) \
				/ np.maximum(np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
				+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
				axis = 1).reshape((-1,1)), 1e-20),0)) - 1/2.) \
				- obs) * cd_vec ),axis= 1) / obs_error ** 2).reshape((-1,1))	
		#------------------------------------------------------------
	else:
		#------------------------------------------------------------
		#calculate residual:
		res = np.sqrt(np.sum(np.square(\
				((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
				+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
				+(px_vec - px_lens_vec) * loc_vec) * c2 * mass * (px_lens_vec - px_vec) \
				/ np.maximum((np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
				+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
				axis = 1).reshape((-1,1)) \
				+ 2 * c2 * mass * (px_lens_vec - px_vec)), 1e-20) \
				- obs) * cd_vec ),axis= 1) / obs_error ** 2).reshape((-1,1))
		#------------------------------------------------------------
	return res 

def FitFunction_Micro_Scandir(par, numtimeloc_vec, sc_scandir, obs, obs_error,cd, exact):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	expect coordinates in mas, and mass in Solar_masses
	par = parameter: mass plus 5 parameter for each star
	numtimeloc_vec = explanation for obs data: Which star at wich epoch from wich position
	obs = observational imput
	cd = cos(dec)
	exact = Bolean, if True use exact formular, else use aproximation

	'''
	c2 = const_Einsteinradius * const_Einsteinradius
	if cd == 0: cd = 1e-20 #avoid devision by 0
	cd_vec = np.array([cd,1])

	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	#order parameter
	mass = par[0]
	astrometry = par[1:].reshape((-1,5))
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

	if exact: 	
		res = np.sum(sc_scandir\
				* (((radec_vec + pm_vec / cd_vec * time_vec  + px_vec * loc_vec) \
				+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
				+ (px_vec - px_lens_vec) * loc_vec)
				* (np.sqrt(np.maximum(1/4. + c2 * mass * (px_lens_vec - px_vec) \
				/ np.maximum(np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
				+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
				axis = 1).reshape((-1,1)), 1e-20),0)) - 1/2.)			
				- obs) * cd_vec),axis= 1) / obs_error
	else:
		res = np.sum(sc_scandir \
				* (((radec_vec + pm_vec / cd_vec * time_vec  + px_vec * loc_vec) \
				+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
				+ (px_vec - px_lens_vec) * loc_vec) \
				* c2 * mass * (px_lens_vec - px_vec) \
				/ np.maximum((np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
				+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
				axis = 1).reshape((-1,1)) \
				+ 2 * c2 * mass * (px_lens_vec - px_vec)), 1e-20) \
				- obs) * cd_vec),axis= 1) / obs_error
	if np.isnan(res.any()): print(1)
	return res 

def Jacobian_Micro(par, numtimeloc_vec, obs, obs_error, cd, exact):
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
	if cd == 0: cd = 1e-20
	cd_vec = np.array([cd,1])
	c2 = const_Einsteinradius * const_Einsteinradius
	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	# order parameters
	mass = par[0]
	astrometry = par[1:].reshape((-1,5))
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

	# calculat partial derivertivs d_delta_ra_i/dxj,d_delta_dec_i/dxj for Motion and Microlensing 
	# cross terms dF(star_i)/ dpar(star_j) for (i != j and j != 0) are equal 0
	Motion = Partial_Motion(time_vec, loc_vec, cd,exact)
	Micro = Partial_Micro(radec_vec, radec_lens_vec, pm_vec, pm_lens_vec,
				px_vec, px_lens_vec, mass, time_vec, loc_vec, cd, exact)
	
	if exact:
		# calculate F_i (see FitFunction)
		f_of_x = np.sqrt(np.sum(np.square(\
			((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
			+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
			+ (px_vec - px_lens_vec) * loc_vec) \
			* (np.sqrt(np.maximum(1/4. + c2 * mass * (px_lens_vec - px_vec) \
			/ np.maximum(np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
			+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
			axis = 1).reshape((-1,1)), 1e-20),0)) - 1/2.) \
			- obs) * cd_vec ),axis= 1) / obs_error **2).reshape((-1,1)) 

		# calculate delta_ra_i,delta_dec_i
		a_d = (((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
			+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
			+ (px_vec - px_lens_vec) * loc_vec) \
			* (np.sqrt(np.maximum(1/4. + c2 * mass * (px_lens_vec - px_vec) \
			/ np.maximum(np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
			+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
			axis = 1).reshape((-1,1)), 1e-20),0)) - 1/2.) \
			- obs) * [cd,1]).reshape((-1,1,2)) / obs_error

	else:	
		# calculate Fi (see FitFunction)

		f_of_x = np.sqrt(np.sum(np.square(\
			((radec_vec + pm_vec / cd_vec * time_vec + px_vec * loc_vec) \
			+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec / cd_vec \
			+(px_vec - px_lens_vec) * loc_vec) * c2 * mass * (px_lens_vec - px_vec) \
			/ np.maximum((np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
			+ (pm_vec - pm_lens_vec) * time_vec + (px_vec - px_lens_vec) * loc_vec * cd_vec), \
			axis = 1).reshape((-1,1)) \
			+ 2 * c2 * mass * (px_lens_vec - px_vec)), 1e-20) \
			- obs) * cd_vec ),axis= 1)).reshape((-1,1)) /obs_error

		# calculate delta_ra_i,delta_dec_i
		a_d = (((radec_vec + pm_vec * time_vec / cd_vec + px_vec * loc_vec) \
			+ ((radec_vec - radec_lens_vec) + (pm_vec-pm_lens_vec) * time_vec /cd_vec \
			+ (px_vec - px_lens_vec) * loc_vec)   \
			* c2 * mass * (px_lens_vec - px_vec) \
			/ np.maximum((np.sum(np.square((radec_vec - radec_lens_vec) * cd_vec \
			+ (pm_vec - pm_lens_vec) * time_vec \
			+ (px_vec - px_lens_vec) * loc_vec * cd_vec), axis = 1).reshape((-1,1)) \
			+ 2 * c2 * mass * (px_lens_vec - px_vec)), 1e-20) \
			- obs) * cd_vec).reshape((-1,1,2)) / obs_error	

	# calculate dFi/dxj
	partial = (1 / np.maximum(f_of_x,1e-20) * (np.sum(a_d * (Motion + Micro) , axis = 2))).reshape(-1)

	# reshape matrix to include all cross terms
	jac = np.zeros(len(par) * len(num_vec))
	ind_par, ind_jac = Jac_Micro_Index(num_vec, len(par))
	jac[ind_jac] = partial[ind_par]

	return jac.reshape(-1, len(par))

def Jacobian_Micro_Scandir(par, numtimeloc_vec, sc_scandir, obs, obs_error, cd, exact):
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


	if cd == 0: cd = 1e-20
	cd_vec = np.array([cd,1])
	c2 = const_Einsteinradius * const_Einsteinradius
	par = np.array(par)
	num_vec = numtimeloc_vec[0]
	time_vec = numtimeloc_vec[1].reshape((-1,1))
	loc_vec = numtimeloc_vec[2]

	# order parameters
	mass = par[0]
	astrometry = par[1:].reshape((-1,5))
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

	# calculat partial derivertivs d_delta_ra_i/dxj,d_delta_dec_i/dxj for Motion and Microlensing 
	# cross terms dF(star_i)/ dpar(star_j) for (i != j and j == 0) are equal 0
	Motion = Partial_Motion(time_vec, loc_vec, cd)

	Micro = Partial_Micro(radec_vec, radec_lens_vec, pm_vec, pm_lens_vec,
						px_vec, px_lens_vec, mass, time_vec, loc_vec, cd, exact)
	# calculate dFi/dxj
	partial = np.sum(sc_scandir.reshape(-1,1,2) * (Motion + Micro) / obs_error.reshape(-1,1,1), axis = 2).reshape(-1) 

	# reshape matrix to include all cross terms

	jac = np.zeros(len(par) * len(num_vec))
	ind_par, ind_jac = Jac_Micro_Index(num_vec, len(par))
	jac[ind_jac] = partial[ind_par]
	if np.isnan(jac.any()): print(2)

	return jac.reshape(-1, len(par))

def Partial_Motion(time_vec, loc_vec, cd):
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
	alpha = np.array([zero, zero, zero, zero, zero, zero, one * cd, zero, time_vec[:,0], zero, loc_vec[:,0] * cd]) 
	delta = np.array([zero, zero, zero, zero, zero, zero, zero, one, zero, time_vec[:,0], loc_vec[:,1]]) 
	return np.array([alpha, delta]).T

def Partial_Micro(radec_vec, radec_lens_vec, pm_vec, pm_lens_vec,
	px_vec, px_lens_vec, mass, time_vec, loc_vec, cd, exact):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	calculat partial derivatives for the Microlensing term
	if exact = True:
		ra_dec = DELTA_ra_dec * [sqrt(1/4+ThetaE**2/Dist**2) - 1/2]

	else:
	 	ra_dec = DELTA_ra_dec*ThetaE**2/(Dist**2+2 * ThetaE**2)

	'''

	c2 = const_Einsteinradius*const_Einsteinradius
	cd_vec = np.array([cd,1])


	#calculate Delta_ra_dec
	dradec = (radec_vec - radec_lens_vec) * cd_vec + (pm_vec - pm_lens_vec) * time_vec  \
			+ (px_vec - px_lens_vec) * loc_vec * cd_vec

	#calculateal Distance 
	Dist2 = dradec[:,0] ** 2 + dradec[:,1] ** 2
	Dist2 = Dist2.reshape(-1,1)
	#calculate Einstein radius
	THETA_E2 = mass * c2 * (px_lens_vec - px_vec)
	THETA_E2 = THETA_E2.reshape(-1,1)
	
	# calculate derivatives 
	if exact:
		Dist2 = np.maximum(Dist2,1e-20)
		SQRT = np.sign(1/4. + THETA_E2 / Dist2)* np.sqrt(np.maximum(np.abs(1/4. + THETA_E2 / Dist2),1e-20))
		#parallel terms
		d_parallel = cd_vec * SQRT - cd_vec / 2 \
			- dradec / SQRT * THETA_E2 / Dist2 ** 2 * dradec * cd_vec
	
		d_pm_parallel = time_vec * SQRT - time_vec / 2 \
			- dradec / SQRT * THETA_E2 / Dist2 ** 2 * dradec * time_vec
	
		#cross terms
		d_cross = -dradec / SQRT * THETA_E2 / Dist2 ** 2  * dradec[:,[0,1]] * cd_vec[::-1] 
		 
		d_pm_cross = -dradec / SQRT * THETA_E2 / Dist2 ** 2 * dradec[:,[0,1]] * time_vec
		
		#parallax
		d_px = loc_vec * SQRT - loc_vec / 2 \
		+ dradec / SQRT * (c2 * mass / Dist2 - 2 * THETA_E2 / Dist2 ** 2 \
		* np.sum(dradec * loc_vec, axis = 1).reshape(-1,1))
	
		#mass
		d_mass = 1/2. *  dradec * c2 * (px_lens_vec-px_vec) / Dist2 / SQRT


	else:
		nominator = np.maximum(Dist2 + 2 * THETA_E2,1e-20)
		
		#parallel terms
		d_parallel = cd_vec * THETA_E2 / nominator \
				-  dradec * dradec * 2 * cd_vec / nominator ** 2 
		
		d_pm_parallel = THETA_E2 * time_vec / nominator \
			-  dradec * dradec * 2 * THETA_E2 * time_vec / nominator ** 2
	
		#cross terms
		d_cross = - 2 * dradec * dradec[:,[0,1]] * cd_vec[::-1] * THETA_E2 / nominator ** 2 

		d_pm_cross = - 2 * dradec * time_vec * dradec[:,[0,1]] * THETA_E2 / nominator ** 2 
		
		#parallax
		d_px =  (loc_vec * THETA_E2 - dradec * c2 * mass) / nominator	\
			- dradec* THETA_E2 * (2 * np.sum(dradec * loc_vec, axis = 1).reshape(-1,1) - 2 * mass * c2) \
			/ nominator ** 2 
		#mass
		d_mass = dradec * Dist2  * c2 * (px_lens_vec-px_vec) / nominator ** 2
	
	#devolve cross and parallel terms

	OneZero = np.array([1,0])
	ZeroOne = np.array([0,1])
	d_ra = d_parallel * OneZero + d_cross * ZeroOne
	d_dec = d_parallel * ZeroOne + d_cross * OneZero
	d_pmra = d_pm_parallel * OneZero + d_pm_cross * ZeroOne
	d_pmdec = d_pm_parallel * ZeroOne + d_pm_cross * OneZero
	

	
	return np.swapaxes(np.array([d_mass, -d_ra, -d_dec, -d_pmra, -d_pmdec, -d_px,
					d_ra, d_dec, d_pmra, d_pmdec, d_px]), 0, 1)

def Jac_Micro_Index(num_vec, len_par):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	#determine which indices in the jacobian are not equal 0 

	i_par= np.arange(len(num_vec)).reshape(-1,1) * 11
	i_jac= np.arange(len(num_vec)).reshape(-1,1) * len_par
	e0 = np.where(num_vec == 0)
	ne0 = np.where(num_vec != 0)
	ind_jac0 = i_jac[e0] + np.arange(1,6) + 5 * num_vec[e0].reshape(-1,1)
	ind_par0 = i_par[e0] + np.arange(6,11)
	ind_jac1 = i_jac[ne0] + np.concatenate((np.arange(0,6) * np.ones((len(ne0[0]),1), int),\
				np.arange(1,6) + 5 * num_vec[ne0].reshape(-1,1)), axis = 1) 
	ind_par1 = i_par[ne0] + np.arange(11)
	ind_jac = np.concatenate((ind_jac0,ind_jac1), axis= None)
	ind_par = np.concatenate((ind_par0,ind_par1), axis= None)
	
	return ind_par,ind_jac
