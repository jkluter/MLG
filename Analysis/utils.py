import numpy as np
import glob

from MLG.Math import percentile
from MLG import path
from MLG.Simulation.MonteCarlo import loadPkl
from MLG.Analysis import Plots
from astropy.table import Table as ATable
from MLG.Analysis import Table

__all__ = ['statistics', 'ALL_Analysis','find_files']

def statistics(MC_Results, relative = False, plot = False, string = 'Event',\
	vary_par = False, v = False, ext_obs = False):
	'''------------------------------------------------------------
		Description: Calculate the persentiles of the distribution 
	---------------------------------------------------------------
		Input:
			MC_Results: Results of the Montecarlo Simulation
			plot: 		Create Histogram of the distribution of values
			string: 	name of the plot's
			vary_par/v: MC_Results contains results for different input parameters
			ext_obs: 	Use the data including external observations  
	---------------------------------------------------------------
		Output:
			stats_all:	List of the statistic results for each event, 
						differen set of input parameters (if vary == TRUE)
						different Type of input parameters    
			Eventnames: List of Names for each event
			Event_out:	Lens, and Source data for each event
	------------------------------------------------------------'''

	#------------------------------------------------------------
	#combine keyword 
	vary_par = any([vary_par,v])
	#------------------------------------------------------------
	
	#------------------------------------------------------------
	#extract data from dictionary
	Events = MC_Results['Eventlist']
	Names = MC_Results['Names']
	if ext_obs:
		Result = MC_Results['Results_ext_obs']
	else: 
		Result = MC_Results['Results']

	# give each event an name if Names is empty
	if Names is None: 
		Names = [str(ev[0].getId()) for ev in Events]
	elif Names == '':
		Names = [str(ev[0].getId()) for ev in Events]
	Names = np.array(Names)

	#remove certain events from outList
	outList = np.array(ATable.read(path + 'out_list', format = 'csv')['source_ID'], str)
	#------------------------------------------------------------
	if relative: 
		#------------------------------------------------------------
		# place holder for doing somthing with relative uncertainties which is not defined
		#------------------------------------------------------------
		pass
	else:
		#------------------------------------------------------------		
		# arrays for storing results, and eventlist cleaning 
		stats_all = []
		Eventnames = []
		Event_out = []
		#------------------------------------------------------------

		#------------------------------------------------------------
		# order the results by determine the number of differnt input parameters
		if vary_par:
			n_par_picks = 0
			m = -1
			for mass in MC_Results['Input_parameter'][0][:,0]:
				if mass != m:
					m = mass
					n_par_picks +=1

		else:
			n_par_picks = 1
		n_per_mass = int(len(Result[0][:,0]) / n_par_picks)
		#------------------------------------------------------------

		# loop over all events 
		for i in range(len(Result)):
			#------------------------------------------------------------
			# check if element in outlist, if so: skip
			if Names[i] in outList: 
				continue
			elif len(Result[i].shape)==1: continue
			elif (Result[i] == -999).all(): continue
			else: Event_out.append(Events[i])

			#------------------------------------------------------------

			#------------------------------------------------------------
			#store results for different sets of inputparameters
			stats_mass = []
			#------------------------------------------------------------

			# loop over input parameters
			for k in range(n_par_picks):
				stats = []
				# loop over parameters
#				if len(Result[i].shape) ==1 : print(i)
				for j in range(Result[i].shape[1]):

					#------------------------------------------------------------
					# determines the type of the parameter
					if j == 0: par_type = 'mass'
					else: par_type = ['px', 'ra', 'dec', 'pmra', 'pmdec'][j % 5]
					if j > 5: par_type = par_type + '_star'+ str(int((j - 1) / 5))
					else: par_type = par_type + '_lens'
					#------------------------------------------------------------

					#------------------------------------------------------------
					# calculate the statistics
					par_result = (percentile(Result[i][(n_per_mass*k): \
							(n_per_mass * (k+1)),j]))
					par_init = (np.mean(MC_Results['Input_parameter'][i][n_per_mass * k: n_per_mass * (k+1),j]))
					accuracy = par_result[-1]/par_init
					#------------------------------------------------------------

					#------------------------------------------------------------
					#create an Histogramm
					if plot:
						#------------------------------------------------------------
						# get Folder
						if len(Events) == 22: 
							Folder = 'Multi_All'
							idd = '' 
						elif len(Events) == 19:	
							idd = '' 
							Folder = 'Multi_Fivepar'	
						elif len(Events) == 14:	
							idd = '' 
							Folder = 'Multi_Sigma'
						elif len(Events) > 100: 
							Folder = 'Single'
							# counter for multiple events of the same lens
							idd = '_%i'%int(np.where(np.where(Names == Names[i])[0] == i)[0])
						else:
							Folder = ''
							# counter for multiple events of the same lens
							idd = '_%i'%int(np.where(np.where(Names == Names[i])[0] == i)[0])
						#------------------------------------------------------------	

						#------------------------------------------------------------	
						# creates Plot 
						if j == 0:
							if vary_par: Namestring = Names[i] +idd+ '_%3f' % par_init
							
							else: Namestring = Names[i] +idd

							Plots.PlotHistSingle(Result[i][n_per_mass * k: \
								n_per_mass * (k+1),j],par_result,par_init, Folder, \
								string + '_' + Namestring)
						#------------------------------------------------------------	

					#------------------------------------------------------------	
					# save statistics for each parameter in a list
					stats.append((i, par_type,par_init,*par_result,accuracy))
				stats_mass.append(stats)
			Eventnames.append(Names[i])
			#------------------------------------------------------------	

			#------------------------------------------------------------
			# stor multiple results for each set of input parameters 
			if vary_par:
				stats_all.append(stats_mass)
			else:
				stats_all.append(stats)
			#------------------------------------------------------------	

	return stats_all, Eventnames, Event_out

def find_files():
	'''------------------------------------------------------------
	Description: Find and sort actual .pkl files 
	---------------------------------------------------------------'''
	ALL_DAT = glob.glob(path+'Data/*.pkl')
	DAT_info = [i.split('/')[-1].split('_') for i in ALL_DAT] #find all *.pkl files
	aa = ['']*12
	for i in DAT_info: 
		which = None
		#-----------------
		# decode File name:
		if i[-3] == '10': #10 years data 
			if i[1]=='single': # single events
				if i[2] =='external': # additional external observation
					if i[0] == 'MC': which = [0,10] # First iteration
					elif i[0] == 'MCGood' : which = [4,11] # vary the input parameters 
				elif i[2] =='vary': # vary_eta 
					continue
				elif i[2] =='test': # vary_eta 
					continue
				elif '-' in i[2]:
					if i[0] == 'MC': which = 0  # First iteration
					elif i[0] == 'MCGood' : which = 4 # vary the input parameters 
			elif i[1] == 'multi': # combined fits of multiple events
				if i[0] == 'MC': which = 1  # First iterartion 
				elif i[0] == 'MCGood': which = 5 #  vary the input parameters 
				else: continue
				if i[2]=='sigma': which +=2 # selection based on sigma_AL
				elif i[2]=='fivepar': which +=1 # selection based on five parameters of the lens	
		elif i[-3] == '5': #5 years data 
			if i[1] == 'single': # single events
				if i[0] == 'MC': # First iteration
					which = 8
				elif i[0] == 'MCGood': # vary the input parameters 
					which = 9
		if which is None: continue
				
		#-----------------
		# save string in the list:
		if not isinstance(which, list):
			which = [which,]
		for which_i in which:
			if aa[which_i] == '': #no element befor
				aa[which_i] = i
			else: #have already an element
				#-----------------
				#check which file  is the most resent one
				tt1 = i[-5].split('-')
				tt2 = aa[which_i][-5].split('-')
				if int(tt1[2]) > int(tt2[2]): aa[which_i] = i # year
				elif int(tt1[2])<  int(tt2[2]): pass
				else:
					if int(tt1[1]) > int(tt2[1]): aa[which_i] = i #month
					elif int(tt1[1])<  int(tt2[1]): pass
					else:	
						if int(tt1[0]) > int(tt2[0]): aa[which_i] = i # date
						else: pass
	return aa

def ALL_Analysis(MC_res10=None, A_res10=None, MC_res5=None, A_res5=None,\
	MC_res_ext=None, A_res_ext=None, onlyTable=False, sort_epoch=True, plot=False):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	aa = find_files()

	#-----------------
	#load files
	print('Load Files')
	if MC_res10 is  None:
		a10 =['_'.join(aa[i]) for i in range(8)]
		MC_res10 = [loadPkl(i) for i in a10]
	if MC_res5 is  None:
		a5 =['_'.join(aa[i]) for i in [8,9]]
		MC_res5 = [loadPkl(i) for i in a5]
	if MC_res_ext is None:
		if aa[10] == aa[0]:
			MC_res_ext = [MC_res10[0],MC_res10[4]]
		else:
			a_ext =['_'.join(aa[i]) for i in [10,11]]
			MC_res_ext = [loadPkl(i) for i in a_ext]


	#-----------------
	# do the analysis
	print('Analysis')
	if A_res5 is  None:
		print('nominal mission')
		A_res5 = [statistics(MC_res5[i], v = i > 0) for i in range(len(MC_res5))]
	if A_res10 is  None:
		print('extended mission')
		if plot: 
			A_res10 = [statistics(MC_res10[i], v = i > 3, plot = i < 4) for i in range(len(MC_res10))]
		else: 
			A_res10 = [statistics(MC_res10[i], v = i > 3, plot = False) for i in range(len(MC_res10))]
	if A_res_ext is  None:
		print('external_obs')
		A_res_ext = [statistics(MC_res_ext[i], v = i > 0, ext_obs = True) for i in range(len(MC_res_ext))]
	if not onlyTable : 
		#-----------------
		print('Create plots:') 	

		print('5years')
		Plots.PlotHistAll(A_res5[0], string = 'mass_accuracy_histogram_5years', extended = False, vary = A_res5[1])
		print('---')
		print('10years')
		Plots.PlotHistAll(A_res10[0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4])
		Plots.PlotHistAll(A_res10[0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4],ty = 'GS',ylim = [0,40])
		Plots.PlotHistAll(A_res10[0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4],ty = 'c',ylim = [0,40])	

		print('---')
		print('future')
		Plots.PlotHistAll(A_res10[0], string = 'mass_accuracy_histogram_future', future = True, vary = A_res10[4]) 
		print('---')
		print('multi_sources')
		Plots.PlotHistAll(A_res10[1],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_all') 
		Plots.PlotHistAll(A_res10[2],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_fps') 
		Plots.PlotHistAll(A_res10[3],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_sig') 		

		Plots.PlotDifferentIdeas(*A_res10[4:])
		Plots.PlotMassAccuracy(A_res10[4])
		#Plots.PlotCHI2(MC_res_ext[0])
		#Plots.PlotMassInOut(A_res10[0])	
	#-----------------
	# Create Table
	print('create_Table single')
	Table.create_Table(A_res10[4],A_res5[1],A_res_ext[1],sort_epoch = sort_epoch)
	print('create_Table multi')
	Table.create_Table_multi(A_res10[4:])
	return (MC_res10, A_res10, MC_res5, A_res5, MC_res_ext, A_res_ext)