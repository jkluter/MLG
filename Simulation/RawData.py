import numpy as np
from MLG.Math import dist, dist_scandir,cosdeg,sindeg,sinmas,cosmas, unitFac
from MLG.StellarMotion import getSun,t_0,stars_position_micro
import time
__all__ = ['Star','Data','sigmaGaia', 'resolvable', 'MassGmag']

class Star(object):
	'''---------------------------------------------------------------
	Defines a star object with 
	position in deg, 
	proper motion in mas/yr
	parallax in mas
	Gmag
	Source_ID
	mass in M_sun
	epoch of the closest aproach in Julian Year
	---------------------------------------------------------------'''
	def __init__(self,pos,pos_err,pm,pm_err,px,px_err, Gmag, mass = 0,mass_err = 0, id = -1, 
		unit_pos = 'deg',unit_pos_err = 'mas', 
		unit_pm = 'mas/yr', unit_px = 'mas',tca = -1, eta = 0):
		'''---------------------------------------------------------------
		initialise the Star
		------------------------------------------------------------------

		Input:
			pos: position vector [ra,dec]
			pos_error: error of the position vector [ra_err,dec_err]
			pm: propermotion vector: [pmra,pmdec]
			pm_err: error of the propermotion vector [ra_err,dec_err]
			px: parallax 
			px_err:  error of the parallax
			Gmag: Magnitude of the Star
			mass: Mass of the lens
			mass_err: error of the mass
			id: source_id
			unit_pos: unit of the given position
			unit_pos_err: unit of the given error of the position
			unit_pm: unit of the given propermotion
			unit_px: unit of the given parallax
			tca:  epoach of the closest approach
			eta: value to store eta 
		---------------------------------------------------------------'''
		self.unit = ['deg','mas/yr','mas','mas']
		self.alpha = pos[0]	* unitFac(unit_pos,'deg')
		self.delta = pos[1] * unitFac(unit_pos,'deg')
		self.alpha_err = pos_err[0] * unitFac(unit_pos_err,'mas')
		self.delta_err = pos_err[1] * unitFac(unit_pos_err,'mas')
		self.pmalpha = pm[0] * unitFac(unit_pm,'mas/yr')
		self.pmdelta = pm[1] * unitFac(unit_pm,'mas/yr')
		self.pmalpha_err = pm_err[0] * unitFac(unit_pm,'mas/yr')
		self.pmdelta_err = pm_err[1] * unitFac(unit_pm,'mas/yr')
		self.px = px * unitFac(unit_pm,'mas')
		self.px_err = px_err * unitFac(unit_pm,'mas')
		self.Gmag = Gmag
		if mass < 0: # use an random value for m
			self.mass = np.random.uniform(0.07,1.5)
			self.mass_err = 0.1*self.mass
		
		else:
			self.mass = mass
			if mass_err == 0 : self.mass_err = 0.1 * mass # use an 10 % error
			else: self.mass_err = mass_err

		self.id = id
		self.tca = tca	
		self.eta = eta

	def getParameters(self, vary_par = 0, **kwargs ):
		'''---------------------------------------------------------------
		returns mass of the lens and the astrometric parameter including variation in the error 
		Input:
			vary_par == 0: return the original values
			vary_par == 1: returns the original values + a pm,px
						from the typical  distribution if pm + px is unknown 
			vary_par == 2: returns the o
			vary_par == 3: returns random mass from the error distrbution + orignal astrometric parameters
			vary_par == 4: returns random mass from the error distrbution + orignal astrometric 
						parameters + random pm,px from the typical distribution if pm + px is unknown 
			vary_par == 5: returns random pics from the error distrbution 
		---------------------------------------------------------------''' 

		#vary the mass
		if vary_par > 2 and self.mass_err > 0: mass = np.random.normal(self.mass, self.mass_err)
		else: mass = self.mass

		#vary all astrometric parameters
		if vary_par % 3 == 2:
			alpha_r =  np.random.normal(self.alpha, self.alpha_err * unitFac(self.unit[3],self.unit[0]))
			delta_r =  np.random.normal(self.delta, self.delta_err * unitFac(self.unit[3],self.unit[0]))
			if self.pmalpha != 0:
				pmalpha_r =  np.random.normal(self.pmalpha, self.pmalpha_err)
				pmdelta_r =  np.random.normal(self.pmdelta, self.pmdelta_err)
				px_r =  np.random.normal(self.px, self.px_err)
			else: 
				#pic from a typical distribution 
				pmtot = np.random.normal(5,3)
				pmdir = np.random.uniform(-np.pi,np.pi)
				pmalpha_r = pmtot * np.cos(pmdir)
				pmdelta_r = pmtot * np.sin(pmdir)
				px_r = np.random.normal(2,1)
				while px_r <= 0:
					px_r = np.random.normal(2,1)
		else:
			alpha_r =  self.alpha
			delta_r =  self.delta
			if self.pmalpha != 0 or vary_par % 3 == 0:
				#orignal astrometric parameters 
				pmalpha_r =  self.pmalpha
				pmdelta_r =  self.pmdelta
				px_r =  self.px
			else: 
				#pic from a typical distribution 
				pmtot = np.random.normal(5,3)
				pmdir = np.random.uniform(-np.pi,np.pi)
				pmalpha_r = pmtot * np.cos(pmdir)
				pmdelta_r = pmtot * np.sin(pmdir)
				px_r = np.random.normal(2,1)
				while px_r <= 0:
					px_r = np.random.normal(2,1)
		return [mass, alpha_r, delta_r, pmalpha_r, pmdelta_r, px_r]

	def getPos(self,unit='deg'):
		#returns the 2015.5 position
		return self.alpha* unitFac('deg',unit),self.delta* unitFac('deg',unit)
	def getPm(self, unit='mas/yr'):
		#returns the propermotion
		return self.pmalpha * unitFac('mas/yr',unit),self.pmdelta * unitFac('mas/yr',unit)
	def getPx(self, unit='mas'):
		#returns the parallax
		try: return self.px * unitFac('mas',unit)
		except:	return self.parallax * unitFac('mas',unit) #old variable name
	def getMag(self):
		#returns the G Magnitude
		return(self.Gmag)
	def getMass(self):
		#returns the Mass of the lens
		return self.mass
	def getId(self):
		#returns the Source_id
		return(self.id)
	def getTca(self): 
		#returns the epoch of the closest approach
		return self.tca
	def getEta(self):
		#returns eta
		return(self.eta)

class Data(object):
	'''---------------------------------------------------------------
	Stars: List of Stars
	NumStars: Number of stars 
	Mass_lens: Mass for each observational Data
	data: [[Star_ID,alpha,delta,Epoch,Scandir, loc, NCCD],...] for each set
	par: input parameters for each set
	par: input parameters for each set
	---------------------------------------------------------------'''
	
	def __init__(self, Stars, Epoch, ScanDir, NCCD, Mass = None,\
		out_unit = 'deg', n_par_picks = None, onlysource = False,timer =False, **kwargs):
		'''---------------------------------------------------------------		
			Creats simulated sets of noice free observational Data (i.e no observational
			errors are included)
			Creats one observation for every source and every epoch (if resolvable)
		------------------------------------------------------------------
			Input:
				Stars: vector of 1 lens and multiple source stars
				Epoch: 	epoch of observation for the event
				ScanDir: Position angle of the scan dirrection for each epoch [rad] 
				NCCD: 	Number of CCD for the observation
				Mass: mass of the lens 
				out_unit: unit of the observations
				n_par_picks2: number of different Sets  
							if none return one set 1 
							if list: [number of differne parameters, number of different masses]
				onlysource: returns only data of the source
				timer: print computing time for different steps)
				**kwargs: vary_par for Star.getParameters 
						  exact for StelarMotion.stars_position_micro
						  ext_obs for resolvable
		---------------------------------------------------------------'''
		#compunting-time-tracker
		cputt = [0,0,0]
		cputt[0] -= time.time()

		#-----------------------------------------------------------------
		# Determine number of different sets
		if n_par_picks is None: 
			n_par = 1
			n_mass = 1
		elif isinstance(n_par_picks, list):
			n_par = n_par_picks[0]	
			n_mass = n_par_picks[1]
		else: 
			n_mass = n_par_picks 
			n_par= n_par_picks
		#-----------------------------------------------------------------
		
		#-----------------------------------------------------------------
		#unit of the data and input parameters 
		self.Stars = np.array(Stars)
		self.unit = [out_unit, *Stars[0].unit]
		# Number of stars 
		self.NumStars = len(Stars)
		# setup Data array
		dat_multi = []  
		# setup array for multiple sets of input parameters
		par_multi = np.zeros((n_par,len(Stars) * 5 + 1))
		# setup array for G magnitudes
		Gmag =np.zeros(len(Stars))
		for pick in range(n_par):
			for i in range(len(Stars)):
				if i == 0: #store parameters of the lens
					par_multi[pick,0:6] = Stars[i].getParameters(**kwargs)
					if Mass is not None: 
						par_multi[pick,0] = Mass
					elif pick%(n_par/n_mass) != 0: 
						par_multi[pick,0] = par_multi[pick-1,0] 
					elif pick == 0:
						self.Mass_lens = par_multi[pick,0] 
				else:
					#store parameters of the Sources
					par_multi[pick,5*i+1 : 5*i+6] = Stars[i].getParameters(**kwargs)[1:]

				# store Gmagnitude of lens and source
				if pick == 0: 
					Gmag[i] = Stars[i].getMag()
		#-----------------------------------------------------------------
		
		#-----------------------------------------------------------------				
		#caluclate sin(alpha),cos(alpha),sin(delta),cos(delta)
		if 'deg'  in self.unit[1]:
			scsc = sindeg(par_multi[0,1]), cosdeg(par_multi[0,1]), sindeg(par_multi[0,2]),\
					max(cosdeg(par_multi[0,2]), 1e-16)
		if 'mas'  in self.unit[1]: 
			scsc = sinmas(par_multi[0,2]), cosmas(par_multi[0,1]), sinmas(par_multi[0,2]),\
					max(cosmas(par_multi[0,2]), 1e-16)
		if 'arcsec'  in self.unit[1]: 
			scsc = sinmas(par_multi[0,1]/1000), cosmas(par_multi[0,1]/1000), sinmas(par_multi[0,2]/1000),\
					max(cosmas(par_multi[0,2]/1000), 1e-16)
		#-----------------------------------------------------------------			
		cputt[0] += time.time()
		cputt[1] -= time.time()

		#-----------------------------------------------------------------
		# positional shift due to the parallax
		# calculat position of Lagrange point 2
		L2 = getSun(t_0 + Epoch).T

		# caluate the invers Jacobean for the parallax effect 
		loc = np.repeat(np.stack([(scsc[0] * L2[:,0] - scsc[1] *  L2[:,1])\
			/ max(scsc[3],1e-6),\
			scsc[1] * scsc[2] *  L2[:,0] + scsc[0] * scsc[2]\
			*  L2[:,1] - scsc[3] *  L2[:,2]]).T.reshape(-1,1,2),self.NumStars,axis = 1)
		#-----------------------------------------------------------------

		#-----------------------------------------------------------------
		# Calculat position for each star
		# epoch for all Stars
		Epoch_vec = np.repeat(Epoch.reshape(-1,1),self.NumStars,axis = 1)
		# identifier for all Stars
		num = np.arange(self.NumStars)
		num_vec = np.repeat(num.reshape(1,-1),len(Epoch), axis = 0)  
		# observed position of lens and sources for every parameterpick
		pos_vec  = stars_position_micro(par_multi,num_vec.reshape(-1),Epoch_vec.reshape(-1),
						loc_vec = loc.reshape(-1,2),scsc = scsc, out_unit = self.unit[0],
						unit_pos = self.unit[1], unit_pm = self.unit[2],
						unit_px = self.unit[3], **kwargs)
		pos_vec = pos_vec.reshape(n_par,-1,self.NumStars,2)
		#-----------------------------------------------------------------
		
		#Translate NCCD and ScanDir into numpy arrays
		ScanDir= np.array(ScanDir)

		NCCD = np.array(NCCD)

		cputt[1] += time.time()
		cputt[2] -= time.time()
		#loop over parameter pics
		for k in range(n_par): 
			#get position for this set of parameters
			pos_k = pos_vec[k]
				
			if onlysource: 
				#return only the observations of the source
				W = np.ones((len(pos_k),len(Gmag)))
				W[:,0]=0
			else:
				#check if lens and sources are resolvable
				which = resolvable(pos_k,Gmag,ScanDir, **kwargs)
			#store ID,alpha,delta,Epoch,ScanDir, Loc, NCCD for each resolved datapoint 
			dat = np.hstack([which[1].reshape(-1,1),pos_k[which],Epoch[which[0]].reshape(-1,1),\
					ScanDir[which[0]].reshape(-1,1),loc[which[0],0,:],NCCD[which[0]].reshape(-1,1)])
			dat_multi.append(dat)
		cputt[2] += time.time()

		#print computing time
		if timer: print('RD: %.2f, %.2f, %.2f' % tuple(cputt))

		self.data= dat_multi
		self.par = par_multi

	def Observations(self, n_error_picks = None, timer = False,**kwargs):
		'''---------------------------------------------------------------
			Description: 
				Creates sets of observations (including noise) 
				from Noise free data
		------------------------------------------------------------------
			Input:
				n_error_picks: number of picks from the errorelipse
		------------------------------------------------------------------
			Output:	
				Obs: list of [alpha,delta] for all Observation
				StarID: list of StarID's for all Observation
				Epoch: list of Epochs for all Observation
				sigAL: list of precision in along-scan dirrection for all Observations
				sc_ScanDir: list of [sin(ScanDir),cos(ScanDir)] for all Observations
		---------------------------------------------------------------'''
		#compunting-time-tracker
		cputt  = [0,0,0,0]
		cputt[0] -= time.time()
	
		#if isinstance(self.data, list):
		#combine all datasets
		dat = np.vstack(self.data)
		#determine length of each dataset 
		length = np.array([len(i) for i in self.data]) 
		# first Index of each dataset 
		index = np.array([np.sum(length[:i]) for i in range(1,len(self.data))]).astype(int)
		#else: 
		#	dat=self.data
	
		if n_error_picks is None:
			n_error_picks = 1
		cputt[0] += time.time()
	
		cputt[1]-= time.time()
		#---------------------------------------------------------------
		#get StarID,Epoch ScanDir NCCD  for each observation
		StarID  = dat[:,0].astype(int)
		Epoch = dat[:,3]
		ScanDir = dat[:,4]
		NCCD = dat[:,7]
		#---------------------------------------------------------------
	
		#---------------------------------------------------------------
		#calculate accuracy of gaia
		if self.Stars is None or len(StarID) == 0:
			#accuracy in along-scan dirrection
			sigAL = sigmaGaia() * np.ones([len(dat[:,0])])
			#set across-scan accuracy to 1 arcsec (is not used for fitting)
			sigAC = 1000 * np.ones([len(dat[:,0])])	
		else:
			#accuracy in along-scan dirrection
			sigAL = sigmaGaia_np(self.Stars[StarID])/np.sqrt(NCCD)
			#set across-scan accuracy to 1 arcsec (is not used for fitting)
			sigAC = 1000 * np.ones([len(dat[:,0])])	
		#---------------------------------------------------------------

		#---------------------------------------------------------------
		#create multidimensional random valu from gaussian distribution
		rand = np.random.normal(0,1,len(dat)*2*n_error_picks).reshape(n_error_picks,-1, 1,2)
		cputt[1] += time.time()
		#---------------------------------------------------------------

		#---------------------------------------------------------------
		#Transform into alpha-delta for each observation
		cputt[2] -= time.time()
		co = np.cos(ScanDir)
		si = np.sin(ScanDir)
		sc_ScanDir = np.vstack([si,co]).T
		cd = cosdeg(dat[:,2])
		trans = np.zeros((1,len(dat),2,2))	
		trans[0,:,0,0] = sigAL * si/cd/3.6e6
		trans[0,:,0,1] = - sigAC * co/cd/3.6e6
		trans[0,:,1,0] = co*sigAL/3.6e6
		trans[0,:,1,1] = si*sigAC/3.6e6  
		Obs = np.sum(trans*rand,3) + dat[:,1:3]
		cputt[2] += time.time()
		#---------------------------------------------------------------

		#---------------------------------------------------------------
		# order output 
		cputt[3] -= time.time()
		#if isinstance(self.data, list):
		Obs = np.split(Obs,index, axis = 1)
		StarID = np.split(StarID,index)
		Epoch =np.split(Epoch,index)
		sigAL = np.split(sigAL,index)
		sc_ScanDir = np.split(sc_ScanDir,index)
		cputt[3] += time.time()
		#---------------------------------------------------------------
		if timer:
			print('OB: %.4f,%.4f,%.4f,%.4f' % tuple(cputt))	
	
		return  Obs,StarID,Epoch, sigAL, sc_ScanDir

def sigmaGaia(Stars=None, Gmag = None):
	'''---------------------------------------------------------------
	Calculate Gaias accuracy in along scan direction from G magnitud
	------------------------------------------------------------------
	Input:
		Stars: Star Object
		Gmag: apparent G magnitude (only if Stars == None)
	------------------------------------------------------------------	
	Output:
		sigmaCCD: Sigma in alongscan direction for one CCD observation in mas  

	---------------------------------------------------------------'''
	if Stars is None:
		if Gmag is not None:
			z = 10 ** (0.4 * (np.maximum(Gmag, 14) - 15))
			sigmaCCD =((-1.631 + 680.766 * z + 32.732 * z**2)**0.5/1000 *7.75+0.1)

		else: return 1 
	else:
			z = 10 ** (0.4 * (np.maximum(Stars.getMag(), 12) - 15))
			sigmaCCD =((-1.631 + 680.766 * z + 32.732 * z**2)**0.5/1000 *7.75+0.1)
	return sigmaCCD
sigmaGaia_np = np.vectorize(sigmaGaia)

def resolvable(pos_vec, Gmag , ScanDir, ext_obs = False,**kwargs):
	'''	---------------------------------------------------------------
	checks if two or more stars are within one readout window 
	Input
		pos_vec: position for each star and each observation 
		Gmag: G Magnitude of the Stars
		ScanDir: Position angle of the scan direction for each observation
		ext_obs: external observations: exclude the last 2*n data points from checking
	------------------------------------------------------------------
	Output
		Indicees for resolvabel observation 
	---------------------------------------------------------------'''
	#-----------------------------------------------------------------
	# limits for the allong-scan and crosscan dirrections direction
	AL_lim= 59*6 *np.ones(len(Gmag))
	AL_lim[np.where(Gmag < 13)] = 59*9 # bright sources hav an larger readout window
	AC_lim= 177*6 *np.ones(len(Gmag))
	Lim = np.vstack([AL_lim,AC_lim]).T
	#-----------------------------------------------------------------
	unitfactor = 3.6e6 #translate deg to mas
	#-----------------------------------------------------------------

	if len(pos_vec.shape) == 3: # check if multiple observations are given
		#split pos to ra and dec
		ra = pos_vec[:,:,0] 
		dec = pos_vec[:,:,1]

		Gmag = np.repeat([Gmag],len(ra),axis = 0) # repeat Gmag for each observation
		
		# [sin cos] and [cos -sin] of the scan dirrection
		sc_scandir = np.array([np.sin(ScanDir),np.cos(ScanDir)]).T 
		cs_scandir = np.array([np.cos(ScanDir),-np.sin(ScanDir)]).T

		
		#Check if all other stars are outside of therir readout window
		Window  = np.sum(np.floor(np.abs(((np.repeat(ra.reshape(len(ra),-1,1,1),len(ra[0]),axis = 2)\
			- np.repeat(ra.reshape(len(ra),1,-1,1),len(ra[0]),axis = 1))\
			*sc_scandir.reshape(-1,1,1,2) * np.cos(dec[:,0]).reshape(-1,1,1,1)\
			+(np.repeat(dec.reshape(len(dec),-1,1,1),len(dec[0]),axis = 2) \
			- np.repeat(dec.reshape(len(dec),1,-1,1),len(dec[0]),axis = 1))\
			*cs_scandir.reshape(-1,1,1,2))*unitfactor/Lim.reshape(1,1,-1,2))),axis =3)
		
		# Check if on star is brighter than the other
		CompG = np.maximum(np.repeat(Gmag.reshape(len(ra),1,-1),len(ra[0]),axis = 1) \
			- np.repeat(Gmag.reshape(len(ra),-1,1),len(ra[0]),axis = 2) -1, 0)
		
		# Identity matrix to avoid comparison with itself
		I = np.repeat([np.eye(len(ra[0]))],len(ra), axis = 0)
		
		for i in range(0,2 * ext_obs):
				I[-i-1] = np.ones(len(ra[0]))
		

		resolve = np.prod(Window+CompG+I, axis = 2) #exclude observation if all are False 
	else: #only one observation is given
		
		#split pos to ra and dec
		ra = pos_vec[:,1]
		dec = pos_vec[:,1]

		# [sin cos] and [cos -sin] of the scan dirrection
		sc_scandir = [np.sin(scandir),np.cos(scandir)]
		cs_scandir = [np.cos(scandir),-np.sin(scandir)]

		#Check if all other stars are outside of therir readout window
		Window = np.sum(np.floor(np.abs(((np.repeat(ra.reshape(-1,1,1),len(ra),axis = 1) \
			- np.repeat(ra.reshape(1,-1,1),len(ra),axis = 0)) * np.cos(dec[0])*sc_scandir\
			+ (np.repeat(dec.reshape(-1,1,1),len(ra),axis = 1 ) \
			- np.repeat(dec.reshape(1,-1,1),len(ra),axis = 0)) * cs_scandir)*unitfactor/Lim)),\
				axis =2)
		
		# Check if on star is brighter than the other
		CompG = np.maximum(np.repeat(Gmag.reshape(1,-1),len(ra),axis = 0) \
				- np.repeat(Gmag.reshape(-1,1),len(ra),axis = 1) -1, 0)
		
		# Identity matrix to avoid comparison with itself
		I = np.eye(len(ra))
		resolve = np.prod(Window+CompG+I, axis = 1)
	return np.where(resolve)

def MassGmag(Gmag,px):
	'''---------------------------------------------------------------
	Returns the mass of an star based on the absolut G magnitude
	and Mass luminosity relations (see Klüter et al 2018)
	---------------------------------------------------------------
	Input 
		Gmag: aparent G magnitude
		px: Parallax  in mas
	Output
		mass: in M_sun
	---------------------------------------------------------------'''
	Gabs = Gmag + 5 * math.log(px/100, 10) 
	a1,b1,c1,a2,b2,n = [7.86232823e-03,  -2.90891912e-01,   1.18248766e+00,  -3.01175062e-01, 1.88952947e+00, 8.71127084e+01]
	Mass =  pow(pow(np.exp(a1*Gabs*Gabs+b1*Gabs +c1),-n)+pow(np.exp(a2*Gabs+b2),-n), -1/n)
	return max(0.07, Mass)