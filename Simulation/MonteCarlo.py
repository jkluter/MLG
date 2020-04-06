'''
This File Provides the Funktions for the Montecarlo simulation
To start use MLG.Start
'''


from MLG import	path,default_table, message_folder,  Version,Subversion
from  MLG.Simulation import RawData
from  MLG.Simulation import RealData
from  MLG.Modeling  import fitting_micro , fitting_motion
from  MLG.Math import percentile
from  joblib import Parallel, delayed
import numpy as np 
import time
import pickle
import datetime
import glob
import os
__all__ = ['StartMC','loadPkl', 'calcParallel']

def StartMC(EventList, namelist = None, keywords = {}, **kwargs):
	'''------------------------------------------------------------
		Description: 
		This Function distribute the calculation for the MonteCarlo
		simulation on multiple CPU_Cores and stors the results in an PKL file	
		uses the joblib package for parallel computing 
	---------------------------------------------------------------
		Input:
			EventList: 	String of the input table 	or 
						List of Lens source paris for each event or
						MC_results for vary the mass also
			namelist: List of names for each event (for ploting routine, if None use of the index)
	---------------------------------------------------------------
		Output:
			MC_result:  Dictnary of the results. Contains:
			 	'Header': List of all control imputs, and Code Version
			 	'Results': List of the fit parameters
			 	'Input_parameter': List of the input parameters
			 	'Results_ext_obs':  List of the fit parameters with external_obs if ext_obs > 0 
			 	'Chi2':	List of the reduced CHI2
			 	'Eventlist': (see Input)
			 	'Names': (see Input)
	------------------------------------------------------------'''
	
	#--------------------------------------------------------------
	#Update controle keywords
	keywords.update(kwargs) 
	num_core = keywords.get('num_core', 6) # Number of cores
	instring = keywords.get('instring', '') # string for the save files
	message = keywords.get('message', False) # write the process of the computation to an file 
	DR3 = keywords.get('DR3', False) # Use of the Data of DR3 only
	extended = keywords.get('extended', False) # Use of the 10 years data of the extended mission
	vary_eta = keywords.get('vary_eta', False) # Use an new scaninglaw for 2020.5-2024.5
	ext_obs		 = keywords.get('ext_obs', False) # include external observations
	n_error_picks = keywords.get('n_error_picks', 500) # Number of pickt Observation from the error ellips
	n_par_picks = keywords.get('n_par_picks', 1) # Number of pickt input Parameters from the error ellips
	namelist = keywords.get('namelist', namelist) # Check if namelist is given in the keyword dictionary 
	#--------------------------------------------------------------


	#--------------------------------------------------------------
	#Create random seeds for the calculation
	seed = []
	for i in range(num_core): 
		seed.append((int(time.time()*10**9))**4 %4294967295)
	keywords['seed']=seed
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	#Load Table if EventList is an string
	if isinstance(EventList, str) == True:
		if EventList == '': 
			EventList,_ = RealData.loadRealData(default_table)
		else:
			EventList,_ = RealData.loadRealData(EventList)
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	# start calculations with varing the mass if a dict is given 
	if isinstance(EventList, dict) == True:
		print('vary mass') 

		#--------------------------------------------------------------
		# extract lists from Dictionary
		MC_Results = EventList 
		res_all_events = MC_Results['Results'] 
		par_all_events = MC_Results['Input_parameter']
		EventList = MC_Results['Eventlist']
		namelist = MC_Results['Names']
		header = MC_Results['Header']
		#--------------------------------------------------------------

		#--------------------------------------------------------------
		# only consider the events with an error < 100%
		EventList = [EventList[i] for i in range(len(res_all_events)) \
				if (percentile(res_all_events[i][:,0])[0]> 0)]
		if namelist is None:
			namelist_good = None
		else:
			namelist_good = [namelist[i] for i in range(len(res_all_events)) \
				if (percentile(res_all_events[i][:,0])[0]> 0)]
		#--------------------------------------------------------------
		
		#--------------------------------------------------------------
		# update control keywords
		keywords['Good']=True # indication calculation of the good events only 
		goodstr = 'Good_' # string for indication calculation of the good events only 
		keywords['vary_par']=5 # vary all parameters 
		n_error_picks = keywords.get('n_error_picks', 500) #check if value is given in keywords else set to defaut 
		keywords['n_par_picks']=n_par_picks 
		n_par_picks = keywords.get('n_par_picks', 100) #check if value is given in keywords else set to defaut
		keywords['n_error_picks']=n_error_picks
		#--------------------------------------------------------------
	
	#--------------------------------------------------------------
	# start first calculation 
	elif isinstance(EventList, list) == True:
		#--------------------------------------------------------------
		# update control keywords
		keywords['n_par_picks']=n_par_picks # Update keyword
		keywords['n_error_picks']=n_error_picks # Update keyword
		keywords['vary_par']=keywords.get('vary_par', 1) # vary only non given parameters 
		keywords['Good']=False # indication calculation of the good events only 
		goodstr=''
		#--------------------------------------------------------------

		#--------------------------------------------------------------
		#set default namelist to integer (not comparable within different inputs)
		if namelist == None: 
			namelist == [str(i) for i in range(len(EventList))]
		#--------------------------------------------------------------
	#--------------------------------------------------------------
	# exclude different filetypes
	else:
		print ('Input Error!')
		return
	#--------------------------------------------------------------
		
	#--------------------------------------------------------------
	# create header 
	if instring is not '':
		instring = instring + '_'
	
	if DR3:
		#use only events before 2019.5 for DR3
		EventList = [EventList[kkk] for kkk in range(len(EventList)) if EventList[kkk][0].getTca() < 2019.5 ]
		header = len(EventList),3,n_error_picks,n_par_picks, keywords,Version+'.'+Subversion
	elif extended or vary_eta or ext_obs:
		#use the data of the extended mission
		header = len(EventList),10,n_error_picks,n_par_picks,keywords, Version+'.'+Subversion
	else:
		header = len(EventList),5,n_error_picks,n_par_picks,keywords, Version+'.'+Subversion
	print(time.ctime())
	print(header)
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	# Distribute on different cores 
	num_events = len(EventList)
	events_per_core = num_events/num_core
	#calculation of multiple events (proxima centauri tend to take as much computation time as all other events)
	if len(EventList[0]) > 10 and num_core > 2: 
		events_per_core = (num_events-1)/(num_core -1)
		for core in range(num_core):
			partstring = path+'Simulation/evpart/eventlist_'+goodstr+instring+'part_%i.pkl' % core 
			if core == 0:
				f = open(partstring, 'wb')
				pickle.dump([EventList[0],], f)
				f.close()
			elif core == num_core -1: 				
				f = open(partstring, 'wb')
				pickle.dump(EventList[1 + round(events_per_core * (core - 1)):], f)
				f.close()
			else: 
				f = open(partstring, 'wb')
				pickle.dump(EventList[1 + round(events_per_core * (core - 1)) : 1 + round(events_per_core * (core))], f)
				f.close()
	#distribute events equaly
	else:
		for core in range(num_core):
			partstring = path+'Simulation/evpart/eventlist_'+goodstr+instring+'part_%i.pkl' % core 
			if core == num_core -1: 				
				f = open(partstring, 'wb')
				pickle.dump(EventList[round(events_per_core * core):], f)
				f.close()
			else: 
				f = open(partstring, 'wb')
				pickle.dump(EventList[round(events_per_core * core) : round(events_per_core * (core + 1))], f)
				f.close()
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	#start calculations parallel
	if num_core != 1:
		res_par = Parallel(n_jobs=num_core)(delayed(calcParallel)(i,instring, keywords) for i in range(num_core))
	else:
		res_par = [calcParallel(0,instring, keywords),]
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	#merge the results for the parallel computations
	res_all_events = []
	par_all_events = []
	res_no_ML_events = []
	chi2_events = []
	for res_par_core in res_par:
		for par_event in res_par_core[1]:
			par_all_events.append(par_event)
		for res_event in res_par_core[0]:
			res_all_events.append(res_event)
		for res_no_event in res_par_core[2]:
			res_no_ML_events.append(res_no_event)
		for res_no_event in res_par_core[3]:
			chi2_events.append(res_no_event)
	if ext_obs:
		MC_Results = {'Header': header,'Results':res_all_events,'Input_parameter':par_all_events,\
			 'Results_ext_obs':res_no_ML_events, 'Chi2':chi2_events,'Eventlist':EventList,'Names':namelist}
	else:
		MC_Results = {'Header': header,'Results':res_all_events,'Input_parameter':par_all_events,\
			 'Results_no_ML':res_no_ML_events, 'Chi2':chi2_events,'Eventlist':EventList,'Names':namelist}
	
	#--------------------------------------------------------------

	#--------------------------------------------------------------
	# save results as pkl file 
	string = path + 'Data/MC'+goodstr[:-1] +'_'+ instring + datetime.datetime.today().strftime('%d-%m-%Y')\
			 + '_%.f_%.f_%.f_%.f.pkl' % (header[:4])
	f = open(string, 'wb')
	pickle.dump(MC_Results, f)
	print(string)
	if message: 
		os.system('cp ' + string + ' ' + message_folder)
	return MC_Results

def loadPkl(filename = '',Good = False,  extended = False, n_error_picks = False, n_par_picks=False):
	'''------------------------------------------------------------ 
		Load the MC_results PKL files ()
	---------------------------------------------------------------
		Input: 
			filename: filename expected in the MLG/Data Folder
	---------------------------------------------------------------
		Output: 
			MC_Results: Dictonary containg results from StartMC
	------------------------------------------------------------'''
	if len(glob.glob(filename)) == 0:
		if Good: good = 'Good' 
		else: good = ''    
		if extended: ex = string(extended) + '_'
		else: ex = '*_'
		if n_error_picks: er = string(extended) + '_'
		else: er = '*_'	
		if n_par_picks: pa = string(extended) + '.pkl'
		else: pa = '*.pkl'
		gstring= (path + 'Data/MC' + good + '*_' + ex + er + pa)
		g = glob.glob(gstring)
		string = g[-1]
		if filename	!= '':
			if len(glob.glob(path + 'Data/' + filename)) == 0:
				print('File not found! Using standard file')
			else: 
				string = glob.glob(path + 'Data/' + filename)[0]
		else:
			print('Using standard file')

	else: string =  glob.glob(filename)[0]
	print(string)
	f = open(string,'rb')
	pkl = pickle.load(f)
	f.close()
	if isinstance(pkl,dict): #from Version 3.1
		if 'Results_comp' in pkl.keys():
			pkl['Results_ext_obs'] = pkl.pop('Results')
			pkl['Results'] = pkl.pop('Results_comp') 
		MC_Results = pkl
	else: 
		#until Version 3.1
		if len(pkl) == 6:
			chi2_events = None
			EventList, res_all_events, par_all_events,res_no_ML_events,namelist,header = pkl
			try:
				par_all_events[0].x
				print(1)
			except AttributeError: pass
			else:
				qq = par_all_events
				par_all_events = res_no_ML_events
				res_no_ML_events = par_all_events
		elif len(pkl) == 7:
			EventList, res_all_events, par_all_events,res_no_ML_events,\
				chi2_events,namelist,header = 	pkl
			try:
				par_all_events[0].x
				print(1)
			except AttributeError: pass
			else:
				qq = par_all_events
				par_all_events = res_no_ML_events
				res_no_ML_events = qq
		# Transform format to 3.1 		
		MC_Results = {'Header': header,'Results':res_all_events,'Input_parameter':par_all_events,\
		 'Results_no_ML':res_no_ML_events, 'Chi2':chi2_events,'Eventlist':EventList,'Names':namelist}
	return MC_Results

def calcParallel(part, instring, keywords = {} ,**kwargs):
	'''------------------------------------------------------------
		Description: 
			creats sets of observations for a part of the Eventlist
			and fits the data
	---------------------------------------------------------------
		Input: 
			part: 	which part of the EventList
			instring: file string of the EventList
			keywords/kwargs: setup values for the simulation (see below)
	---------------------------------------------------------------
		Output:
			res_part: results of the individual fits (contains a separate list for each events)
			par_part: Inputparameters of the indiciual fits (contains a separate list for each events)
			res_no_ML: result without microlensing 
			chi2_red_part: Chi2 values of the individual fits (contains a separate list for each events)
	------------------------------------------------------------'''
	
	#---------------------------------------------------------------
	#extract keywords
	keywords.update(kwargs)
	seed		 = keywords.get('seed', None) # seed's for the randomised process
	Good		 = keywords.get('Good', False) # Second loop with variation of the mass

	extended	 = keywords.get('extended', False) # Use of the data for extended Mission 
	ext_obs		 = keywords.get('ext_obs', False) # ext_obs include external observations
	DR3		 	 = keywords.get('DR3', False) # Use of the data for DR3 only 
	exact		 = keywords.get('exact', False)	# Use the exact astrometric shift or the approximation

	n_error_picks= keywords.get('n_error_picks', 500) # number of picks from the error elips of the measurements
	n_par_picks  = keywords.get('n_par_picks', 1) # number of different picks of input parameters
	vary_par 	 = keywords.get('vary_par', 1) # which parameters should be varied for n_par_picks
	
	prefit_motion= keywords.get('prefit_motion', False) # fit the propermotion with out microlensing as preput
	vary_eta 	 = keywords.get('vary_eta', False) #Use data while vary eta
	onlysource	 = keywords.get('onlysource', False) # use the data of the source only

	timer		 = keywords.get('timer', False) # print computing-time for different steps
	message		 = keywords.get('message', False) # save messeage file for keeping track of the process 
	silent		 = keywords.get('silent', True) # boolen if information shoul not be printed
	#---------------------------------------------------------------

	#---------------------------------------------------------------
	# computationtime tracker
	cputt = [0.,0.,0.,0.,0.,0.,0]
	cputt[0] -= time.time()
	# update seed for the randomised process
	if seed is not None:
		np.random.seed(seed[part])
	#---------------------------------------------------------------
	
	#---------------------------------------------------------------
	# initilise arrays for storing the results 
	# fit results
	res_part = []
	# fit results without microlensing
	res_no_ML = []
	# input parameters
	par_part = []
	# Chi2_reduced value 
	chi2_red_part = []
	chi2_red_mot_part = []
	#---------------------------------------------------------------

	#---------------------------------------------------------------
	#	load part of the Eventlist
	if Good: 
		partstring = path+'Simulation/evpart/eventlist_Good_'+instring+'part_%i.pkl' % part 	
		f = open(partstring,'rb')
		EventList = pickle.load(f)
		f.close
		os.remove(path+'Simulation/evpart/eventlist_Good_'+instring+'part_%i.pkl' % part)

	else:	
		partstring = path+'Simulation/evpart/eventlist_'+instring+'part_%i.pkl' % part 
		f = open(partstring,'rb')
		EventList = pickle.load(f)
		f.close
		os.remove(path+'Simulation/evpart/eventlist_'+instring+'part_%i.pkl' % part)
	cputt[0] += time.time()
	#---------------------------------------------------------------

	for i1 in range(len(EventList)):
		#---------------------------------------------------------------
		# updat process tracker
		if message: 
			os.system('touch %s.%s-%s-0.message'%(message_folder+instring[:-1],part,i1))
		#---------------------------------------------------------------

		cputt[1] -= time.time()

		#---------------------------------------------------------------
		# initilise list for each event 
		# fit results
		res_single = []
		res_external_obs = [] #for comparison (test)
		# input parameters
		par_single = []
		# Observations (for fitting  without microlensing)
		Obs_save = []
		# Chi2_reduced value 
		chi2_red_single = []
		#---------------------------------------------------------------

		#get Stars
		Stars = EventList[i1]
		
		#---------------------------------------------------------------
		# create observations from real scanning law
		RealObs = RealData.getRealObs(Stars[0],Stars[1],**keywords) 
		Data = RawData.Data(Stars,*RealObs,**keywords)
		Obs,Starnumber,Epoch,Obs_err, sc_scandir = Data.Observations(**keywords)
		#---------------------------------------------------------------
		
		cputt[1] += time.time()

		for i2 in range(n_par_picks): 
			#loop over each set of differen input parameters
			par = Data.par[i2]
			for i3 in range(n_error_picks):
				#loop over each set of differen picks from the error ellips of the Observation

				cputt[6] = -time.time()
				cputt[2] -= time.time()
				

				# include outliers 5% chance 
				if ext_obs:
					outliers = np.random.uniform(0,1,len(Starnumber[i2]))
					outliers[-4*ext_obs:] = 1 #keep the external observations
					n = np.where(outliers > 0.05)
				else:
					n = np.where(np.random.uniform(0,1,len(Starnumber[i2])) > 0.05)

				#---------------------------------------------------------------
				# get the right set of observations
				Obs_i = Obs[i2][i3][n]
				if len(Obs_i) < 11+4*ext_obs:
					# not enought observation 
					re = -999* np.ones(len(Stars)*5+1)
					res_single.append(re)
					chi2_red_single.append(-999.)
					par_single.append(par)
					cputt[2] += time.time()
					res_external_obs.append(re)

				else:
					# get the right set of observations
					# (identical within each set of input parameters)
					Starnumber_i = Starnumber[i2][n]
					Epoch_i = Epoch[i2][n]
					Obs_err_i = Obs_err[i2][n]
					sc_scandir_i = sc_scandir[i2][n]
					#---------------------------------------------------------------
					
					#---------------------------------------------------------------
					#compare with external observations (test)
					if ext_obs:
						res_ext= fitting_micro.Fit_Micro(Obs_i,Starnumber_i,Epoch_i,Obs_err_i, sc_scandir_i,\
						 		bounds = True, **keywords)
						res_external_obs.append(res_ext.x)

						# remove external observations
						Obs_i=Obs_i[:-4*ext_obs]
						Obs_err_i = Obs_err_i[:-4*ext_obs]
						Starnumber_i = Starnumber_i[:-4*ext_obs]
						Epoch_i = Epoch_i[:-4*ext_obs]
						sc_scandir_i = sc_scandir_i[:-4*ext_obs]

					#---------------------------------------------------------------


					#---------------------------------------------------------------
					# store observations
					Obs_save.append([Obs_i, Starnumber_i,Epoch_i, Obs_err_i, sc_scandir_i])
					cputt[2] += time.time()
					#---------------------------------------------------------------

					#---------------------------------------------------------------
					# Fit model to the observational data
					cputt[3] -= time.time()
					res= fitting_micro.Fit_Micro(Obs_i,Starnumber_i,Epoch_i,Obs_err_i, sc_scandir_i,\
						 bounds = True, **keywords)
					cputt[3] -= time.time()
					#---------------------------------------------------------------

					#---------------------------------------------------------------
					#store results
					cputt[4] -= time.time()
					if len(res.x) != len(Stars)*5+1:
						#Check if all stars have observations
						re = -999* np.ones(len(Stars)*5+1)
						re[:len(res.x)] = res.x
						re[len(res.x):] = None
						if len(res.x) == 6 : re[0] = -999
						res_single.append(re)
					else:
						res_single.append(res.x)	
					chi2_red_single.append(res.cost * 2 / (len(res.fun) - 11))
					par_single.append(par)
					cputt[4] += time.time()
					#---------------------------------------------------------------



				#---------------------------------------------------------------
				#print computing time (for optimisation)
				if timer:
					if not (i3 % 100): 
						cputt[5] += time.time()
						print('Part %i-%i Step %i %.2f' % (part,i1,i3,cputt[6]))
						if time.time() - ti > 10: return 0 
				#---------------------------------------------------------------
			
			#---------------------------------------------------------------
			# updat process tracker
			if message: 
				os.system('mv %s.%s-%s-%s.message %s.%s-%s-%s.message'%\
						(message_folder + instring[:-1], part, i1, i2,\
						message_folder + instring[:-1], part, i1, i2+1))
			#---------------------------------------------------------------


		if not silent: print(i3, time.time() - ti)
		#---------------------------------------------------------------
		# sort results
		cputt[5] -= time.time()
		res_single = np.stack(res_single)
		res_part.append(res_single)
		par_single = np.stack(par_single)
		par_part.append(par_single)
		chi2_red_single = np.stack(chi2_red_single)
		chi2_red_part.append(chi2_red_single)
		#---------------------------------------------------------------
		#compare without external observations (test)
		if ext_obs:
			res_external_obs = np.stack(res_external_obs)
			res_no_ML.append(res_external_obs)

		#---------------------------------------------------------------

		#---------------------------------------------------------------
		# fit mean result also without microlensing 
		if len(res_single[:,0]) > 1:
			if ext_obs:
				chi2_red_micro_best = 0
				chi2_red_motion = 0
			else:	
				p50 = percentile(res_single[:,0])[1]
				if p50 == -999:
					chi2_red_mot_part.append(-999)
					chi2_red_micro_best = -999
					res_no_ML.append(-999*np.ones(5*len(Stars)))
					chi2_red_motion = -999
				else:
					best = np.where(res_single[:,0] == p50)[0][0]
					res_motion = fitting_motion.Fit_Motion(*Obs_save[best],**keywords)
					res_no_ML.append(res_motion.x)
					chi2_red_micro_best = chi2_red_single[best] 
					chi2_red_motion= res_motion.cost * 2 / (len(Obs_save[best][1]) - 10)
					chi2_red_mot_part.append(chi2_red_motion)


			cputt[5] += time.time()
		#---------------------------------------------------------------

		#---------------------------------------------------------------
		# print computing time (for optimisation)
		if timer: 
			print('MC %.0f-%.0f: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f' % (part,i1,*tuple(cputt)))
		#---------------------------------------------------------------

		else:
			#---------------------------------------------------------------
			# print statistics 
			if len(res_single[:,0]) > 1: 
				print('%.0f-%.0f: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %s'\
					% (part,i1, Stars[0].mass,*percentile(res_single[:,0]), chi2_red_micro_best,chi2_red_motion, time.ctime().split(' ')[3]))
			#---------------------------------------------------------------

	if not silent: print('Part {} done!'.format(part))
	
	#---------------------------------------------------------------
	# updat process tracker
	if message:
		os.system('rm %s.%s*.message'%(message_folder+instring[:-1],part))
		os.system('touch %s.%s.Done.message'%(message_folder+instring[:-1],part))
	#---------------------------------------------------------------

	return res_part, par_part, res_no_ML,chi2_red_part
