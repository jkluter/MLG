import sys
import time

'''------------------------------------------------------------
	Description: 
---------------------------------------------------------------
	Input:
---------------------------------------------------------------
	Output:
------------------------------------------------------------'''

__all__ = ['start']



def start(ev,setup_keys):
	import MLG
	from MLG.Simulation.MonteCarlo import StartMC

	keywords = {}
	keywords['exact']= False
	keywords['num_core'] = setup_keys.get('num_core', 6) # Number of cores
	instring = setup_keys.get('instring', None) # string for the save files
	keywords['message'] = setup_keys.get('message', False) # write the process of the computation to an file 
	keywords['extended'] = setup_keys.get('extended', False) # Use of the 10 years data of the extended mission
	vary = setup_keys.get('vary_mass', True) 
	if 'ext_obs' in setup_keys.keys():
		ext_obs = True
		keywords['ext_obs'] = setup_keys.get('ext_obs', 0)
	else:
		ext_obs = False

	keywords['n_error_picks'] = setup_keys.get('n_error_picks', 500) # Number of pickt Observation from the error ellips
	keywords['n_par_picks']  = 1 # Number of pickt input Parameters from the error ellips
	if 'DR3' in setup_keys.keys(): 
		keywords['DR3'] = setup_keys.get('DR3')
	if 'vary_eta' in setup_keys.keys(): 
		keywords['vary_eta'] = setup_keys.get('vary_eta')
	if instring is None:
		if setup_keys['type'] >= 9: instring = ''
		else: instring = ['','single','multi_all','multi_fivepar',\
			'multi_all, multi_fivepar','multi_sigma',\
			'multi_all, multi_sigma','multi_fivepar, multi_sigma',\
			'multi_all, multi_fivepar, multi_sigma', ][setup_keys['type']]
		if ext_obs: instring = instring + '_external_obs' 
		elif keywords['extended']  == False:   instring = instring+ '_5years'
		if setup_keys.get('best',False): instring=instring + '_best'
	keywords['instring']  = instring
	if isinstance(ev,dict):
		print('dict')
		res1 = ev
	else:
		ti = time.time()
		res1 = StartMC(ev, keywords=keywords)
		print(time.time() -ti)
	if vary:
		keywords['n_par_picks'] = setup_keys.get('n_par_picks', 100) # Number of pickt input Parameters from the error ellips
		ti = time.time()
		res2 = StartMC(res1, keywords=keywords)
		print (time.time() - ti)	



class GUI_Ana:
	color_hide = '#808080'
	color_vis = '#ffffff'

	def __init__(self, window12):
		aa = MLG.Analysis.find_files()
		self.Analysis = False
		self.keywords = {}	
		self.window12 = window12
		self.window12.title('Gui')
		self.window12.geometry('800x700')
		self.MC_res = None
		self.Analysis_res = None

		#------------------------------------------------------------
		# Single
		self.line()
		tk.Label(self.window12, text="fixed parameters", borderwidth=2, relief="flat")
		tk.Label(self.window12, text="fixed parameters", borderwidth=3, relief="flat")
		tk.Label(self.window12, text="Single 10yrs:").grid(row = 11, column = 0, sticky=W)
		self.MC_res10_all = tk.Entry(self.window12)
		self.MC_res10_all.insert(END, '_'.join(aa[0]))
		self.MC_res10_all.grid(row=11,column=1, sticky=W,columnspan=8)

		tk.Label(self.window12, text="Single 5yrs:").grid(row = 12, column = 0, sticky=W)
		self.MC_res5_all = tk.Entry(self.window12)
		self.MC_res5_all.insert(END, '_'.join(aa[8]))
		self.MC_res5_all.grid(row=12,column=1, sticky=W,columnspan=8)

		tk.Label(self.window12, text="Single ext:").grid(row = 13, column = 0, sticky=W)
		self.MC_res10_ext = tk.Entry(self.window12)
		self.MC_res10_ext.insert(END, '_'.join(aa[10]))
		self.MC_res10_ext.grid(row=13,column=1, sticky=W,columnspan=8)

		# single _ vary 
		self.MC_res10_all_v = tk.Entry(self.window12)
		self.MC_res10_all_v.insert(END, '_'.join(aa[4]))
		self.MC_res10_all_v.grid(row=11,column=3, sticky=W,columnspan=4)

		self.MC_res5_all_v = tk.Entry(self.window12)
		self.MC_res5_all_v.insert(END, '_'.join(aa[9]))
		self.MC_res5_all_v.grid(row=12,column=3, sticky=W,columnspan=4)

		self.MC_res10_ext_v = tk.Entry(self.window12)
		self.MC_res10_ext_v.insert(END, '_'.join(aa[11]))
		self.MC_res10_ext_v.grid(row=13,column = 3, sticky=W,columnspan=4)
		#------------------------------------------------------------

		#------------------------------------------------------------
		# Multiple
		self.line()
		tk.Label(self.window12, text="fixed parameters", borderwidth=2, relief="flat")
		tk.Label(self.window12, text="fixed parameters", borderwidth=3, relief="flat")
		tk.Label(self.window12, text="Multi all:").grid(row = 21, column = 0, sticky=W)
		self.Multi_all = tk.Entry(self.window12)
		self.Multi_all.insert(END, '_'.join(aa[1]))
		self.Multi_all.grid(row=21,column=1, sticky=W,columnspan=8)

		tk.Label(self.window12, text="Multi 5par:").grid(row = 22, column = 0, sticky=W)
		self.Multi_5par = tk.Entry(self.window12)
		self.Multi_5par.insert(END, '_'.join(aa[2]))
		self.Multi_5par.grid(row=22,column=1, sticky=W,columnspan=8)

		tk.Label(self.window12, text="Multi sig:").grid(row = 23, column = 0, sticky=W)
		self.Multi_sig = tk.Entry(self.window12)
		self.Multi_sig.insert(END, '_'.join(aa[3]))
		self.Multi_sig.grid(row=23,column=1, sticky=W,columnspan=8)

		# Multiple_vary 
		self.Multi_all_v = tk.Entry(self.window12)
		self.Multi_all_v.insert(END, '_'.join(aa[5]))
		self.Multi_all_v.grid(row=21,column=3, sticky=W,columnspan=4)

		self.Multi_5par_v = tk.Entry(self.window12)
		self.Multi_5par_v.insert(END, '_'.join(aa[6]))
		self.Multi_5par_v.grid(row=22,column=3, sticky=W,columnspan=4)

		self.Multi_sig_v = tk.Entry(self.window12)
		self.Multi_sig_v.insert(END, '_'.join(aa[7]))
		self.Multi_sig_v.grid(row=23,column = 3, sticky=W,columnspan=4)
		self.line()
		#------------------------------------------------------------

		self.button1 = tk.Button(self.window12, text='Load_Data', command = self.Load,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=30, column = 0,columnspan = 1)	

		self.button1 = tk.Button(self.window12, text='Analysis', command = self.Do_Analysis,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=30, column = 1,columnspan = 1)	

		self.button1 = tk.Button(self.window12, text='Create Plot', command = self.Plot_All,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=30, column = 2,columnspan = 1)	

		self.button1 = tk.Button(self.window12, text='Create Table', command = self.Create_Table,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=30, column = 3,columnspan = 1)	

		self.PlotHistAll = tk.IntVar()
		self.PlotDifferentIdeas = tk.IntVar()
		self.PlotMassAccuracy = tk.IntVar()
		self.PlotCHI2 = tk.IntVar()
		self.PlotMassInOut = tk.IntVar()
		tk.Label(self.window12, text="Histogram Mass").grid(row = 40, column = 0, sticky=W)
		tk.Label(self.window12, text="Multi Different Cases ").grid(row = 40, column = 1, sticky=W)
		tk.Label(self.window12, text="Sigma M vs M").grid(row = 40, column = 2, sticky=W)
		tk.Label(self.window12, text="Histogram Chi^2").grid(row = 40, column = 3, sticky=W)
		tk.Label(self.window12, text="Mass_out vs Mass_in").grid(row = 40, column = 4, sticky=W)
		tk.Checkbutton(self.window12, variable=self.PlotHistAll ).grid(row=41, column = 0, sticky=W)
		tk.Checkbutton(self.window12, variable=self.PlotDifferentIdeas ).grid(row=41, column = 1, sticky=W)
		tk.Checkbutton(self.window12, variable=self.PlotMassAccuracy).grid(row=41, column = 2, sticky=W)
		tk.Checkbutton(self.window12, variable=self.PlotCHI2).grid(row=41, column = 3, sticky=W)
		tk.Checkbutton(self.window12, variable=self.PlotMassInOut).grid(row=41, column = 4, sticky=W)

		self.line()	
		self.button1 = tk.Button(self.window12, text='        Exit        ', command = self.Exit,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=200, column = 1,columnspan = 2)	

	def line(self):
		tk.Label(self.window12, text="-------------------------------------------------------------------------------------------------------").grid(column =0 ,columnspan=4)	
	def Load(self):
		print('Load_Data')
		files_10 = [self.MC_res10_all.get(),self.Multi_all.get(),self.Multi_5par.get(),self.Multi_sig.get(),self. MC_res10_all_v.get(),self.Multi_all_v.get(),self.Multi_5par_v.get(),self.Multi_sig_v.get()]
		files_5 = [self.MC_res5_all.get(),self.MC_res5_all_v.get()]
		files_ext =[self. MC_res10_ext.get(),self.MC_res10_ext_v.get()]
		MC_res10 = [MLG.Simulation.MonteCarlo.loadPkl(i) for i in files_10]
		MC_res5 = [MLG.Simulation.MonteCarlo.loadPkl(i)  for i in files_10]
		MC_ext = [MLG.Simulation.MonteCarlo.loadPkl(i)  for i in files_ext]
		self.MC_res = [MC_res10,MC_res5,MC_ext]
	def Do_Analysis(self):
		if self.MC_res is None: self.Load()
		print('Analysis')
		self.Analysis_res = []
		for j in range(3):
			if j ==0: print('10years')
			if j ==1: print('5years')
			if j ==2: print('external_obs') 
			self.Analysis_res.append([MLG.Analysis.statistics(self.MC_res[j][i], v = i >= len(self.MC_res[j])/2, ext_obs = (j==2)) for i in range(len(self.MC_res[j]))])

	def Plot_All(self):
		if self.Analysis_res is None: self.Do_Analysis()
		print('Create plots:') 	
		print(self.PlotHistAll.get())
		if self.PlotHistAll.get():
			print('5years')
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[1][0], string = 'mass_accuracy_histogram_5years', extended = False, vary = A_res5[1])
			print('---')
			print('10years')
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4])
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4],ty = 'GS',ylim = [0,40])
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][0], string = 'mass_accuracy_histogram_10years', vary = A_res10[4],ty = 'c',ylim = [0,40])	
			print('---')
			print('future')
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][0], string = 'mass_accuracy_histogram_future', future = True, vary = A_res10[4]) 
			print('---')
			print('multi_sources')
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][1],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_all') 
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][2],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_fps') 
			MLG.Analysis.Plots.PlotHistAll(self.Analysis_res[0][3],xlim = [0,1.3], string = 'mass_accuracy_histogram_multi_sig') 		
		if self.PlotDifferentIdeas.get(): 
			string = MLG.Analysis.Plots.PlotDifferentIdeas(*self.Analysis_res[0][4:])
			for i in strings:
				root = tk.Tk()  
				canvas = tk.Canvas(root, width = 500, height = 500)  
				canvas.pack()  
				image = tk.Image.open(string)  
				image = image.resize((480, 440))
				img = tk.ImageTk.PhotoImage(image)
				canvas.create_image(20, 20, anchor=NW, image=img) 
				root.mainloop() 
		if self.PlotMassAccuracy.get(): 
			MLG.Analysis.Plots.PlotMassAccuracy(self.Analysis_res[0][4])
		if self.PlotCHI2.get(): 
			MLG.Analysis.Plots.PlotCHI2(self.MC_res[0][0])
		if self.PlotMassInOut.get(): 
			MLG.Analysis.Plots.PlotMassInOut(self.Analysis_res[0][0])	

	def Create_Table(self):
		if self.Analysis_res is None: self.Do_Analysis()
		MLG.Analysis.Table.create_Table(self.Analysis_res[0][4],self.Analysis_res[1][1],self.Analysis_res[2][1],sort_epoch = sort_epoch)
		MLG.Analysis.Table.create_Table_multi(self.Analysis_res[0][4:])


	def Exit(self):
		self.window12.quit()
	
	

	def test_values(self):
		self.nCore.delete(0,END)
		self.nCore.insert(0, '4')
		self.nPar.delete(0,END)
		self.nPar.insert(0, '2')
		self.nError.delete(0,END)
		self.nError.insert(0, '10')	
	def check(self):
		if self.StartOK.get() ==  True:
			return True,self.Analysis, self.ev,self.keywords
		else:
			return False,self.Analysis,False,False

	
	#if StartOK.get(): 
	#	start(ev['ev'],keywords)	


class GUI:
	color_hide = '#808080'
	color_vis = '#ffffff'
	default_stringes = ['single', 'single_external_obs', 'single_5years','multi_all','multi_fivepar','multi_sigma', '']	

	def __init__(self, window):
		self.Analysis = False
		self.keywords = {}	
		self.window = window
		self.window.title('Gui')
		self.window.geometry('650x700')
		self.ev_type = tk.IntVar()
		self.ev_type.set(1)
		self.multi_all = tk.IntVar() 
		self.multi_5par = tk.IntVar() 
		self.multi_sig = tk.IntVar() 	

		self.extobs = tk.IntVar()
		self.extended = tk.IntVar()
		self.extended.set(1)
		self.approx= tk.IntVar()
		self.approx.set(0)
		self.vary_mass = tk.IntVar()
		self.message = tk.IntVar()
		self.test = tk.IntVar()

		self.StartOK = tk.IntVar()

		self.line()
		self.L1 = tk.Label(window, text="List of events", borderwidth=2, relief="flat")
		self.L1.grid(row=10, column =0 ,columnspan=4)   	

		tk.Label(self.window, text="Single events:").grid(row = 11, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.ev_type,onvalue=1,command=self.ev_type_c).grid(row=11, column = 1, sticky=W)
		tk.Label(self.window, text="Multiple events:").grid(row = 12, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.ev_type,onvalue=2,command=self.ev_type_c).grid(row=12, column = 1, sticky=W)	
		self.lall = tk.Label(self.window, text="All Sources:")
		self.call = tk.Checkbutton(self.window, variable=self.multi_all,command=self.ev_type_c)
		self.l5par = tk.Label(self.window, text="Sources with pm:")
		self.c5par = tk.Checkbutton(self.window, variable=self.multi_5par,command=self.ev_type_c)
		self.lsig  = tk.Label(self.window, text="Sources with sigma < 0.5:")
		self.csig  = tk.Checkbutton(self.window, variable=self.multi_sig,command=self.ev_type_c)
		self.call.grid(row=12, column = 3, sticky=W)
		self.c5par.grid(row=13, column = 3, sticky=W)
		self.csig.grid(row=14, column = 3, sticky=W)
		self.lall.grid(row = 12, column = 2, sticky=W)
		self.l5par.grid(row = 13, column = 2, sticky=W)
		self.lsig.grid(row = 14, column = 2, sticky=W)	

		tk.Label(self.window, text="Other:").grid(row = 15, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.ev_type,onvalue=3,command=self.ev_type_c).grid(row=15, column = 1, sticky=W)	

		self.otherEV = tk.Entry(self.window)
		self.otherEV.insert(END, 'EventTable.fits')
		self.otherEV.grid(row=15,column=2, sticky=W,columnspan=4)
		self.line()
		self.L2 = tk.Label(self.window, text="File string", borderwidth=0, relief="flat")
		self.L2.grid(row=20, column =0, sticky=W)   
		self.instring = tk.Entry(self.window)
		self.instring.insert(0, self.default_stringes[0])
		self.instring.grid(row=20,column=2, sticky=W,columnspan=4)
		self.line()
		self.L3 = tk.Label(self.window, text="Observations", borderwidth=2, relief="flat")
		self.L3.grid(row=30, column =0 ,columnspan=4)   
		tk.Label(self.window, text="Use extended mission:").grid(row = 31, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.extended,command=self.ev_type_c).grid(row=31, column = 1,columnspan=2, sticky=W)
		tk.Label(self.window, text="Use external observations:").grid(row = 32, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.extobs,command=self.ev_type_c).grid(row=32, column = 1,columnspan=2, sticky=W)
		self.Nextobs = tk.Entry(self.window)
		self.Nextobs.insert(0, '2')
		self.Nextobs.grid(row=32,column=2, sticky=W)
		self.line()
		self.L4 = tk.Label(self.window, text="Monte Carlo Simulation ", borderwidth=2, relief="flat")
		self.L4.grid(row=40, column =0 ,columnspan=4)   	
		
	

		tk.Label(self.window, text="Picks from error elipse:").grid(row = 41, column = 0, sticky=W)
		self.nError = tk.Entry(self.window)
		self.nError.insert(0, '500')
		self.nError.grid(row=41,column=2, sticky=W,columnspan=1)	

		tk.Label(self.window, text="Vary input Parameters:").grid(row = 42, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.vary_mass).grid(row=42, column = 1,columnspan=2, sticky=W)
		self.nPar = tk.Entry(self.window)
		self.nPar.insert(0, '100')
		self.nPar.grid(row=42,column=2, sticky=W,columnspan=1)
		tk.Label(self.window, text="Use center of light :").grid(row = 43, column = 0, sticky=W)
		tk.Checkbutton(self.window, variable=self.approx, onvalue = 0, offvalue = 1).grid(row=43, column = 1,columnspan=2, sticky=W)	

		self.line()
		self.L5 = tk.Label(self.window, text="Parallel computing", borderwidth=2, relief="flat")
		self.L5.grid(row=50, column =0 ,columnspan=4)   	

		tk.Label(self.window, text="Number of Threads ").grid(row = 51, column = 0, sticky=W)
		self.nCore = tk.Entry(self.window)
		self.nCore.insert(0, '16')
		self.nCore.grid(row=51,column=2, sticky=W,columnspan=1)	

		self.line()	
		self.L_error = tk.Label(self.window)

	#button = Button(window, text = 'Start',command = ex)
		self.button1 = tk.Button(self.window, text='        Start        ', command = self.execute,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=100, column = 1,columnspan = 2)	
		self.button2 = tk.Button(self.window, text='Print command', command = self.printCommand,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button2.grid(row=100, column = 3,columnspan = 1)	

		tk.Checkbutton(self.window, text = 'store updates to file', variable=self.message).grid(row=100, column = 0,columnspan=2)	
		tk.Checkbutton(self.window, text = 'test', variable=self.test, command = self.test_values).grid(row=101, column = 0,columnspan=2)

		self.line()	
		self.button1 = tk.Button(self.window, text='        Analysis        ', command = self.start_Analysis,height =2 ,activeforeground = '#0055ff',relief="groove")
		self.button1.grid(row=200, column = 1,columnspan = 2)	


	def line(self):
		tk.Label(self.window, text="-------------------------------------------------------------------------------------------------------").grid(column =0 ,columnspan=4)	
	def ev_type_c(self):
			a = self.instring.get().split(', ')
			if a[0] in self.default_stringes:
				if self.ev_type.get() == 1:
					self.instring.delete(0,END)
					if self.extobs.get():
						self.instring.insert(0, self.default_stringes[1])
					elif self.extended.get(): 
						self.instring.insert(0, self.default_stringes[0])
					else: 
						self.instring.insert(0, self.default_stringes[2])
				if self.ev_type.get() == 2:
					self.instring.delete(0,END)
					c = 0	
					c2= 2
					for nn in self.multi_all.get(), self.multi_5par.get(), self.multi_sig.get():
						c2+=1
						if nn: 
							if c != 0: 
								self.instring.insert(END, ', ')
							c += 1  
							self.instring.insert(END, self.default_stringes[c2])
					if c == 0:
						self.multi_all.set(1)
						self.instring.insert(END, self.default_stringes[3])	

	def start_Analysis(self):
		self.Analysis = True
		self.StartOK.set(0)
		self.window.quit()
	def execute(self):
		import MLG
		from MLG import default_table, inputpath	
		from MLG.Simulation.MonteCarlo import StartMC

		self.L_error.grid_forget()
		ee = self.ev_type.get()
		if ee == 3: 
			if otherEV.get() == 'best' :
				self.keywords['type'] = 1
				eventlist, _ =  MLG.Simulation.RealData.loadRealData(default_table[:-5]+'_best.fits')

			elif '.pkl' in otherEV.get(): 
				try: 
					ev = MLG.Simulation.MonteCarlo.loadPkl(otherEV.get())
				except:
					L_error['text'] ="Can not read " + otherEV.get()
					L_error.grid(row = 99, column = 0, columnspan = 4, sticky=W)
					return
			else:
				self.keywords['type'] = 9
				try: 
					if otherEV.get() == 'abc': pass
					else:
						eventlist, _ =  MLG.Simulation.RealData.loadRealData(otherEV.get(), format=otherEV.get().split('.')[-1])
				except:
					L_error['text'] ="Can not read " + otherEV.get()
					L_error.grid(row = 99, column = 0, columnspan = 4, sticky=W)
					return
			# other
		elif ee == 2:
			# multi
			cc =  self.multi_all.get()+ self.multi_5par.get()*2 + self.multi_sig.get()*4
			self.keywords['type'] = 1 + cc
			strr = [inputpath+str(i)+'.fits' for i in range(1,23)]
			ev_multi = []
			if self.multi_all.get():
				eventlist, _ = MLG.Simulation.RealData.loadRealDataMulti(strr,fivepar = False)
				ev_multi.append(eventlist) 
			if self.multi_5par.get():
				eventlist, _ = MLG.Simulation.RealData.loadRealDataMulti(strr)
				ev_multi.append(eventlist) 
			if self.multi_sig.get():
				eventlist, _ = MLG.Simulation.RealData.loadRealDataMulti(strr,sig = 0.5)
				ev_multi.append(eventlist) 
			if len(ev_multi)>= 1: 
				eventlist = ev_multi
				print('do multiple thinks') 
		elif ee == 1:
			self.keywords['type'] = 1
			eventlist, _ =  MLG.Simulation.RealData.loadRealData(default_table)	
			# single
		else:
			L_error['text'] ="Can not read " + otherEV.get()
			L_error.grid(row = 99, column = 0, columnspan = 4, sticky=W)
			return	

		print(len(eventlist))
		self.keywords['instring'] = self.instring.get()
		#ev = eventlist
		if self.extobs.get(): self.keywords['ext_obs'] = int(self.Nextobs.get())
		else: self.keywords['ext_obs'] = 0
		self.keywords['type']
		self.keywords['extended'] = self.extended.get()
		self.keywords['num_core'] = int(self.nCore.get())
		self.keywords['n_error_picks'] = int(self.nError.get())
		self.keywords['n_par_picks'] = int(self.nPar.get())
		self.keywords['exact']= self.approx.get()
		self.window.quit()
		self.ev = eventlist
		self.StartOK.set(1)
		self.keywords['vary_mass'] = self.vary_mass.get()
	
	def printCommand(self):
		bt = 0
		command_sting = 'python ' + '/'.join(__file__.split('/')[:-1])
		if self.ev_type.get() == 3:
			if self.otherEV.get() == 'best' :
				bt = 1
				tp = 1
			else: tp = 9
		elif self.ev_type.get() == 2:			
			tp = 1+self.multi_all.get()+ self.multi_5par.get()*2 + self.multi_sig.get()*4
		elif self.ev_type.get() == 1: 
			tp =1
		command_sting += ' -j %i'%(tp)
		if tp == 9: command_sting += ' -f ' + self.otherEV.get()
		command_sting += ' -c %s'%(self.nCore.get())
		a = self.instring.get().split(', ') 
		if a[0] in self.default_stringes:
			pass 
		else: 
			command_sting += ' -s %s'%(' '.join(a))

		if int(self.nError.get()) != 500: command_sting +=' -ne %s'%(self.nError.get())
		if self.vary_mass.get(): 
			if int(self.nPar.get()) != 100: command_sting +=' -np %s'%(self.nPar.get())
		else:
			command_sting +=' -v 0'

		if self.extended.get():
			command_sting +=' -e'
		if self.extobs.get():	 command_sting+=' -obs'

		if self.approx.get(): command_sting+=' -exact'
		if bt: command_sting+=' -b'
		if self.message.get(): command_sting+=' -m'
		self.command = tk.Entry(self.window, width=70)
		self.command.insert(0, command_sting)
		self.command.grid(row=103,column=0, sticky=W,columnspan=4)	


	def test_values(self):
		self.nCore.delete(0,END)
		self.nCore.insert(0, '4')
		self.nPar.delete(0,END)
		self.nPar.insert(0, '2')
		self.nError.delete(0,END)
		self.nError.insert(0, '10')	
	def check(self):
		if self.StartOK.get() ==  True:
			return True,self.Analysis, self.ev,self.keywords
		else:
			return False,self.Analysis,False,False

	
	#if StartOK.get(): 
	#	start(ev['ev'],keywords)	

def main():
	error = False
	file = ''
	setup_keys = {}
	for i in range(len(arg)):
		if i == 0: 
			setup_keys['type'] = 1
			setup_keys['ncore'] =6
		elif '-j' in arg[i]: #job_id
			try:
				setup_keys['type'] = int(arg[i+1])
			except ValueError: 
				error = 1
			except IndexError: 
				error = 1
		elif '-f' in  arg[i]:
			try:
				file = arg[i+1]
			except ValueError: 
				error = 2
			except IndexError: 
				error = 2
		elif '-c' in arg[i] or '-ncore' in arg[i]: #number of threads
			try:
				setup_keys['num_core'] = int(arg[i+1])
			except ValueError: 
				error = 3
			except IndexError: 
				error = 3
		elif '-s' in arg[i]:
			string = ''
			for k in range(i,len(arg[i])):
				if arg[k][0] == '-' : break
				else:  	
					if string == '':
						string = arg[k]
					else: string = string + ', '+ arg[k]
			else: setup_keys['instring'] = string  
		elif '-ne' in arg[i]:
			try:
				setup_keys['n_error_picks'] = int(arg[i+1])
			except ValueError: 
				error = 4
			except IndexError: 
				error = 4
		elif '-np' in arg[i]:
			try:
				setup_keys['n_par_picks'] = int(arg[i+1])
			except ValueError: 
				error = 5
			except IndexError: 
				error = 5
		elif '-obs' in arg[i]:
			try:
				setup_keys['ext_obs'] = int(arg[i+1])
			except ValueError: 
				error = 6
			except IndexError: 
				error = 6
		elif '-exact' in arg[i]:	
			setup_keys['exact'] = True
		elif '-e' in arg[i]:	
			setup_keys['extended'] = True
		elif '-DR3' in arg[i]:	
			setup_keys['DR3'] = True
		elif '-vary_eta' in arg[i]:
			setup_keys['vary_eta'] = True
		elif '-v' in arg[i]:	
			try:
				if '-' in arg[i+1]:
					setup_keys['vary_mass'] = True
				else:
					setup_keys['vary_mass'] = (int(arg[i+1]) >= 1)
			except ValueError: 
				error = 7
			except IndexError: 
				setup_keys['vary_mass'] = True		
		elif '-m' in arg[i]:	
			setup_keys['message'] = True
		elif '-b' in arg[i]:	
			setup_keys['best'] = True
	if setup_keys['type'] < 1 or setup_keys['ncore']  < 1 :
		error = 8 
	if error:  
		print('do not understand input %i'%(error))
	else:
		#load MLG libary 
		from MLG import default_table, inputpath, path
		from MLG.Simulation.MonteCarlo import StartMC,loadPkl
		from MLG.Simulation.RealData import loadRealDataMulti,loadRealData
		#load needed input Tables 
		strr = [inputpath+str(i)+'.fits' for i in range(1,23)]
		if '.pkl' in file:
				eventlist = loadPkl(file)
		elif setup_keys['type'] == 1:
			if setup_keys.get('vary_eta', False): eventlist, _ =  loadRealData(inputpath + 'events_eta_3.fits', vary_eta_step = 0.5)
			elif setup_keys.get('best', False): eventlist, _ =  loadRealData(default_table[:-5]+'_best.fits')
			else : eventlist, _ =  loadRealData(default_table)
		elif setup_keys['type']	< 9:
			ev = []
			if (setup_keys['type']-1)%2:
				eventlist, _ = loadRealDataMulti(strr,fivepar = False)
				ev.append(eventlist) 
			if (setup_keys['type']-1)//2%2:
				eventlist, _ = loadRealDataMulti(strr)
				ev.append(eventlist) 
			if (setup_keys['type']-1)//4%2:
				eventlist, _ = loadRealDataMulti(strr,sig = 0.5)
				ev.append(eventlist) 
			if len(ev)>= 1: eventlist = ev
		else: 
			eventlist = loadRealData(file)
		start(eventlist,setup_keys)	

if __name__ == '__main__':
	'''	check non default inputs'''
	arg = sys.argv
	if len(arg) == 1: 
		from tkinter import END,W
		import tkinter as tk
		import tkinter.font as font 
		window = tk.Tk()
		cc = GUI(window)
		window.mainloop()
		GUI_output = cc.check()
		if GUI_output[0] == True:
			start(*GUI_output[2:])
		elif GUI_output[1] == True:
			window.destroy()
			import MLG
			window2 = tk.Tk()
			dd = GUI_Ana(window2)
			window2.mainloop()
	else:
		main()		
else:
	#load MLG libary 
	from MLG import default_table, inputpath
	from MLG.Simulation.MonteCarlo import StartMC
	from MLG.Simulation.RealData import LoadRealDataMulti,LoadRealData

	#load all input Tables 
	eventlist, eventtab =  loadRealData(default_table)
	strr = [inputpath+str(i)+'.fits' for i in range(1,23)]
	eventlist_multi1, eventtab_multi1 = loadRealDataMulti(strr,fivepar = False)
	eventlist_multi2, eventtab_multi2 = loadRealDataMulti(strr)
	eventlist_multi3, eventtab_multi3 = loadRealDataMulti(strr,sig = 1)
	names = [str(i) for i in range(len(eventlist))]


