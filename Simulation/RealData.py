import numpy as np
from astropy.table import Table, vstack
from MLG.Simulation.RawData import Star, sigmaGaia
import MLG.Math 

from MLG import path
from astropy.table import Table,vstack
__all__ = ['getRealObs','varyEtaTab','loadRealData','loadRealDataMulti', 'RealStar']

#---------------------------------------------------------------
#Default table of Gaia observations from GOST
ObsTab  =  Table.read(path+'InputTable/Observations.csv',format = 'csv') 
ObsTabExtended  =  Table.read(path+'InputTable/Observations_extended.csv',format = 'csv') 
ObsTabDR3 = ObsTab[np.where(ObsTab['ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]'] < 2457894)] 

#---------------------------------------------------------------

#---------------------------------------------------------------
#Table of expected observations for differentet eta
ObsTab_vary_eta = []
def varyEtaTab():
	global ObsTab_vary_eta
	#ObsTab_vary_eta1 = Table.read(path+'InputTable/tab_for_vary_eta.csv',format = 'csv')
	ObsTab_vary_eta = Table.read(path+'InputTable/tab3_for_vary_eta.csv',format = 'csv') 
	#ObsTab_vary_eta = vstack([ObsTab_vary_eta1,ObsTab_vary_eta2])
#---------------------------------------------------------------

def getRealObs(star,source = None, obs_string = None, extended = True, DR3 = False, ext_obs = False,vary_eta = False, 
	colnames = ['Target','ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]','scanAngle[rad]','CcdRow[1-7]'], **kwargs):
	'''------------------------------------------------------------
		load real gaia observation form CSV file 
		the CSV file must contain the Source_id of the lens,
		the Epoch of the observation in Julian Date
		the Position angle of the scan direction in rad
		and the number of ccd  
	---------------------------------------------------------------
		Input:
			star: Rawdata.Star objekt for the lensing star 
			source: Rawdata.Star objekt for the source star, used for the resolution of external observation
			obs_string (string): Name of the CSV file. If None use the default table  
			extended (boolean): use of the 10 years mission data
			DR3 (boolean): use the 3 years mission data of Data Releas 3
			ext_obs: include n external observtion with a 7d period
			colnames [string,string,string,string]: List of the columnnames in the CSV table for
												[source_id, T, Scan_direction, CCD_row] 
	
	---------------------------------------------------------------
		Output
			T [Astropy column]:		   of the epoch observations in years after 2015.5
			ScanDir [Astropy column]:  scan dirrection for each observation
			NCCD [Astropy column]:    number of CCD's in the corresponding row
	------------------------------------------------------------'''

	#---------------------------------------------------------------
	#Load Data:
	if obs_string is None: 
		if vary_eta: 
			if len(ObsTab_vary_eta) == 0: varyEtaTab() # load Table scaninglaw.main() results
			eta = star.getEta()//0.05*0.05			# round eta to step size
			Tmin = 2459397.875 						# (2020.5) starting point of the eta variation
			Tab_eta = ObsTab_vary_eta[np.where(np.round(ObsTab_vary_eta['eta[degree]']*20) == round(eta*20))] # select given eta
			if len(Tab_eta) == 0: print(eta)
			Tab = ObsTabExtended[np.where(ObsTabExtended['ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]'] < Tmin)] 
		elif DR3:
			Tab = ObsTabDR3
		elif extended:
			Tab = ObsTabExtended
		else:
			Tab = ObsTab
	else:  Tab = Table.read(obs_string, format ='csv') 
	#---------------------------------------------------------------
	
	#---------------------------------------------------------------
	# creats columns
	# find observations for the given star
	which_obs = np.where(Tab[colnames[0]] == star.getId())

	# get epochs
	T = (Tab[which_obs][colnames[1]] - 2457206.375)/365.25 # Translate Julian Data Years after J2015.5
	
	# translate CCD row to number of CCD's
	NCCD = Tab[which_obs][colnames[3]]
	NCCD[np.where(NCCD != 4)] = 9
	NCCD[np.where(NCCD == 4)] = 8			#Midle row has only 8 CCD's in the astrometric field
	
	# get the scan direction  
	ScanDir = Tab[which_obs][colnames[2]] 	
	#---------------------------------------------------------------

	#---------------------------------------------------------------
	# add 2 perpendicular observations at the epoch of the closest aproach (or a distance 0f 150 mas)
	for dt in range(ext_obs):
		
		# set scan dirrection and number of CCD's
		ScanDir = ScanDir.insert(len(ScanDir),0.)
		ScanDir = ScanDir.insert(len(ScanDir),np.pi/2.)
		NCCD = NCCD.insert(len(NCCD),9)
		NCCD = NCCD.insert(len(NCCD),9)
		
		# Epoch of the closest approach
		tobs = star.getTca()-2015.5

		if source is not None and False:	
			# calculate distance between lens and source as function of time 
			radec_lens =np.array(star.getPos()) 
			pmradec_lens = np.array(star.getPm()) * np.array([1/MLG.Math.cosdeg(radec_lens[1]),1])
			radec_source =np.array(source.getPos()) 
			pmradec_source = np.array(source.getPm()) * np.array([1/MLG.Math.cosdeg(radec_lens[1]),1])
			dist = lambda t: (MLG.Math.dist(radec_lens + t*pmradec_lens/3.6e6, radec_source + t*pmradec_source/3.6e6, unit = 'mas'))
				
			if dt%2:
				while dist(tobs) < 150: 	#limits the impact parameter to 150 mas
					tobs = tobs + 1/365.25
			else:
				while dist(tobs) < 150: 	#limits the impact parameter to 150 mas
					tobs = tobs - 1/365.25
		T = T.insert(len(T),tobs + (dt+1)//2 * (dt%2-0.5)*2/365.25)
		T = T.insert(len(T),tobs + (dt+1)//2 * (dt%2-0.5)*2/365.25)
	#---------------------------------------------------------------


	#---------------------------------------------------------------
	# combine the observations befor 2020.5 with observation for different values of eta
	if  vary_eta:
		# find observations for the given star
		which_obs_eta = np.where(Tab_eta[colnames[0]] == star.getId())
		# get epochs 
		T_eta = (Tab_eta[which_obs_eta][colnames[1]] - 2457206.375)/365.25 

		# translate CCD row to number of CCD's
		NCCD_eta = Tab_eta[which_obs_eta][colnames[3]]
		NCCD_eta[np.where(NCCD_eta != 4)] = 9
		NCCD_eta[np.where(NCCD_eta == 4)] = 8			#Midle row has only 8 CCD's in the astrometric field
		
		# get the scan direction  
		ScanDir_eta = Tab_eta[which_obs_eta][colnames[2]] 	

		#combine observations 
		T=vstack([T,T_eta]) [colnames[1]]
		NCCD=vstack([NCCD,NCCD_eta]) [colnames[3]]
		ScanDir=vstack([ScanDir,ScanDir_eta])[colnames[2]]
	#---------------------------------------------------------------
	return T,ScanDir, NCCD


def loadRealDataMulti(table,format='fits', fivepar = True, sig  = 10):
	'''---------------------------------------------------------------
	load real data for events with multiple background stars from a plc subtables (only on lens with multiple background sources in each table)
	------------------------------------------------------------------
	Input
		table [AstropyTable,...] or [string,...]: loaded tables or paths for the plc subtables
		format [string]:  format of the table if a string is given
		fivepar (boolean): use events with a 5-parameter solution of the background source.
		sig (float): use events with a accuracy of the individual measurements of the background source better than sig 
	------------------------------------------------------------------
	Output
		Stars [list]: list of lens and source pairs
		Tab[AstropyTable]: loaded Table
	---------------------------------------------------------------'''

	#Check if only one table is given
	if not isinstance(table,list): table = [table,]

	#Storage for Starlist and tables 
	Stars_multi = []
	Tab_multi = []

	for i in range(len(table)):
		#Load plc Table 
		if isinstance(table[i],str):
			Tab = Table.read(table[i],format=format) 
		else:
			Tab = table[i]

		Stars= [RealStar(Tab[0],lens=True),] #Create RawData.Star objekt for the lens
		
		#go thrue the list of backgroud sources
		for TabRow in Tab:
			if fivepar:
				Source  = RealStar(TabRow) #Create RawData.Star objekt for the source
				if Source.px_err != 0 and sigmaGaia(Source)/3 < sig: 
					Stars.append(Source) #Add to list of  sources 
			else:	
				Source  = RealStar(TabRow) #Create RawData.Star objekt for the source
				if sigmaGaia(Source)/3 < sig :
					Stars.append(Source) #Add to list of  sources 

		if len(Stars) > 1 : #Check if there is at least one background source
			Stars_multi.append(Stars) #Add to output list
			Tab_multi.append(Tab) #Add to output list

	return Stars_multi,Tab_multi

def loadRealData(table,format='fits',vary_eta_step = 0):
	'''---------------------------------------------------------------
	load real data from a plc subtables (one event per row)
	------------------------------------------------------------------
	Input
		table [AstropyTable or string]: loaded table or path of plc subtables
		format [string]:  format of the table if a string is given
		vary_eta_step [float]: step size of varying eta, if != 0 
	------------------------------------------------------------------
	Output
		Stars [list]: list of lens and source pairs
		Tab[AstropyTable]: loaded Table
	---------------------------------------------------------------'''
	
	#Load plc Table if not already loaded
	if isinstance(table,str):
		Tab = Table.read(table,format=format)
	else:Tab = table

	#---------------------------------------------------------------
	#loop over Tab lines
	Stars = []
	for i in range(len(Tab)):
		if vary_eta_step > 0:
			for eta in np.arange(0, 360, vary_eta_step):
				Stars.append([RealStar(Tab[i], lens=True, eta=eta), RealStar(Tab[i])]) #Create RawData.Star objekt for the lens and source
		else:
			Stars.append([RealStar(Tab[i],lens=True),RealStar(Tab[i])]) #Create RawData.Star objekt for the lens and source

	return Stars,Tab

def RealStar(tabrow,lens=False, eta = 0):
	'''---------------------------------------------------------------
	Creates a RawData.Star object from a given table row
	------------------------------------------------------------------
	Input:
		tabrow:  Row of an Plc Table contains the source_id, pos,pm, px, gmag, mass and epoch data 
		lens: 	select data of the lens
		eta: 	value for eta if vary_eta  	
	------------------------------------------------------------------
	Output 
		Star:   RawData.Star
	---------------------------------------------------------------'''
	
	
	#Check input type
	if 'ra' in tabrow.colnames: #Gaia like table with 'ob_' for source
		
		if lens:LensSource = '' #get lens Data
		else:   LensSource = 'ob_' #get source Data

		#---------------------------------------------------------------
		# get values

		# position
		alpha = tabrow[LensSource + 'ra']
		alpha_err = tabrow[LensSource + 'ra_error']
		delta = tabrow[LensSource + 'dec']
		delta_err = tabrow[LensSource + 'dec_error']
		unit_pos = 'deg'
		unit_pos_err = 'mas'
		# proper motion
		pmalpha = tabrow[LensSource + 'pmra']
		pmalpha_err = tabrow[LensSource + 'pmra_error']
		pmdelta = tabrow[LensSource + 'pmdec']
		pmdelta_err = tabrow[LensSource + 'pmdec_error']
		unit_pm = 'mas/yr'
		# parallax
		parallax = tabrow[LensSource + 'parallax']
		parallax_err = tabrow[LensSource + 'parallax_error']
		unit_px = 'mas'

		#check all 5 parameter solution  exist if not set values to 0
		if pmalpha == 0 or np.isnan(pmalpha):
			pmalpha = 0
			pmalpha_err = 0
			pmdelta = 0
			pmdelta_err = 0
			parallax = 0
			parallax_err = 0

		#G Magnitude
		Gmag = tabrow[LensSource + 'phot_g_mean_mag']
		
		#mass of the Lens
		if lens:
			mass = tabrow['mass']
			mass_error = tabrow['mass_error']
		else: 
			mass = 0
			mass_error = 0

		#epoch of the closest approach
		tca = tabrow['date']
		#---------------------------------------------------------------

	else: #original PLC data (GAVO) 'lens_' for lens 'ob_' for source
		if lens:LensSource = 'lens_'
		else:   LensSource = 'ob_'
		
		#---------------------------------------------------------------
		# get values
		#position
		alpha = tabrow[LensSource + 'ra']
		alpha_err = tabrow[LensSource + 'err_ra']
		delta = tabrow[LensSource + 'dec']
		delta_err = tabrow[LensSource + 'err_dec']
		unit_pos = 'deg'
		unit_pos_err = 'deg'

		# proper motion
		pmalpha = tabrow[LensSource + 'pmra']
		pmalpha_err = tabrow[LensSource + 'err_pmra']
		pmdelta = tabrow[LensSource + 'pmdec']
		pmdelta_err = tabrow[LensSource + 'err_pmdec']
		unit_pm = 'deg/yr'

		# parallax
		parallax = tabrow[LensSource + 'parallax']
		parallax_err = tabrow[LensSource + 'err_parallax']
		unit_px = 'deg'

		#check all 5 parameter solution  exist if not set values to 0
		if pmalpha == 0 or np.isnan(pmalpha):
			pmalpha = 0
			pmalpha_err = 0
			pmdelta = 0
			pmdelta_err = 0
			parallax = 0
			parallax_err = 0
		#G Magnitude
		Gmag = tabrow[LensSource + 'phot_g_mean_mag']
		
		#Mass of the lens
		if lens:
			mass = tabrow['mass']
			mass_error = tabrow['err_mass']
		else: 
			mass = 0
			mass_error = 0

		#epoch of the closest approach
		tca = tabrow['tca']
		#---------------------------------------------------------------
	#---------------------------------------------------------------

	#---------------------------------------------------------------
	#return a Star object
	return Star([alpha, delta], [alpha_err, delta_err], [pmalpha, pmdelta], [pmalpha_err, pmdelta_err],\
				parallax, parallax_err, Gmag, mass, mass_error, tca=tca, id=tabrow[LensSource + 'source_id'],\
				unit_pos= unit_pos, unit_pos_err=unit_pos_err, unit_pm=unit_pm, unit_px=unit_px, eta=eta)
	#---------------------------------------------------------------

