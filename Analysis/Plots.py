
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from MLG.utils import color_own 
from MLG.name import name, name_good 
from MLG.Math import percentile, cosdeg, dist
from MLG import imagepath, paperpath
from MLG.Microlensing.const import const_Einsteinradius

from  MLG.Simulation import RealData
from astropy.table import Table
import os
import time

__all__ = ['dark','PlotHistSingle','PlotDifferentIdeas','PlotMassAccuracy',\
		'PlotMassInOut', 'PlotHistAll','PlotMinMout','PlotCHI2']

def dark(onof = 0):
	if onof is 'on': plt.style.use('dark_background')
	elif onof is 'off': plt.style.use('default')
	elif onof == True: plt.style.use('dark_background')
	else: plt.style.use('default')

def PlotHistSingle(dat,dat_stat,dat_init,Folder= '', string = 'Event', par_type = 'Mass'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
		Plot histogram for an given fit parameter
	'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
		color =   color_own([0.,1,1,1])
	else: 
		black = color_own([0.,0.,0.,1])
		color =   color_own([0.,1,1,1])
	rng = dat_stat[1] + 4 * dat_stat[3], dat_stat[1] + 4 * dat_stat[4]
	fig = plt.figure(figsize = (12,10))
	#fig = plt.figure(figsize = (7,6))

	ax = fig.add_subplot(111)
	plt.subplots_adjust(
		top=0.9,
		bottom=0.15,
		left=0.13,
		right=0.92,
		hspace=0.2,
		wspace=0.2)	

	k = plt.hist(dat,bins = 25,range = rng, color = color, rwidth=0.9)
	y = [max(k[0]) / 2,]*3
	plt.plot(dat_stat[0:3], y, 'x-',color = color_own([1,0.5,0,1]), linewidth = 5, mew = 6, markersize = 30)
	plt.plot((dat_init, dat_init), (0, max(k[0])), '-',color= color_own([1,0,0,1]), linewidth = 10)
	u = plt.ylim()
	plt.ylim([u[0], 1.1* u[1]])
	if par_type == 'Mass':
		plt.text(dat_init,u[1]*1.05,'$M = %.2f\,M_{\odot}$'%dat_init, verticalalignment = 'top',\
					horizontalalignment = 'center', fontsize = 40)



	if par_type == 'Mass': plt.xlabel('$M\,[M_{\odot}]$', fontsize = 50)
	plt.ylabel('#', fontsize = 50)
	ax.tick_params(axis = 'both', which = 'major', labelsize = 30)
	ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
	plt.xticks( fontsize = 40)
	plt.yticks( fontsize = 40)

	# Save files
	if len(Folder) != 0: 
		if Folder[-1] != '/': Folder =  Folder + '/'

	if not os.path.isdir(imagepath +Folder): os.mkdir(imagepath+ Folder)

	if not os.path.isdir(imagepath +Folder+ 'Good_Events/'): 
			os.mkdir(imagepath + Folder + 'Good_Events/')
	if not os.path.isdir(imagepath +Folder+'All_Events/'): 
			os.mkdir(imagepath + Folder +'All_Events/')
	if dat_stat[0] > 0:
		#fig.savefig(imagepath + 'Good_Events/' + string)
		for acc in [15,30,50,100]:
			if (dat_init-dat_stat[0])/dat_init < acc/100.:
				axrange = ax.get_xlim()
				center = (axrange[0]+axrange[1])/2
				#fac = 2 * (100-acc)/85+ 1 * (acc-15)/85 
				if acc == 15:
					fac = 3
				elif acc == 30:
					fac = 2
				elif acc == 50:
					fac = 1.5

				else: fac=1 

				axrange_new = [(axrange[0]-center)*fac + center, \
							(axrange[1] - center) * fac + center]
				plt.xlim(axrange_new)
			


				if not os.path.isdir(imagepath + Folder+'Good_Events/p%i/'%acc):
							os.mkdir(imagepath + Folder+'Good_Events/p%i/'%acc)
				fig.savefig(imagepath + Folder+'Good_Events/p%i/'%acc + string)

				break
	if string == 'Event_2390377345808152832_0':
			plt.xlim([-1,2])
			fig.savefig(paperpath + 'Event_10')
	if string == 'Event_4293318823182081408_0':
			fig.savefig(paperpath + 'Event_70')
	if string == 'Event_4451575895403432064_0':	
			fig.savefig(paperpath + 'Event_30')	
	if string == 'Event_470826482635701376_3':	
			fig.savefig(paperpath + 'Event_400')
	fig.savefig(imagepath + Folder+ 'All_Events/' + string)

	plt.close(fig)

def PlotDifferentIdeas(dat_single,dat_multi_all,dat_multi_fps,dat_multi_sig, \
	string = 'Diff_ID', k = None, Folder = 'Multi_events/'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
		cc = color_own([[0.6,1.2,0,0.9],[0,1,1,0.9],[1,1,0,0.9] ,[1,0,1,0.9]])
	else: 
		black = color_own([0.,0.,0.,1])
		cc = color_own([[0,1,0,0.9],[0,1,1,0.9],[1,.5,0,0.9] ,[1,0,0,0.9]])
	t = [0,0,0,0,0,0]

	#match_events betwenn differnt methods
	ev = []
	c=0
	try :
		a = dat_single[0][0][0][0][0]
	except:
		error = 1
	else: 
		error = 3 
	ti = time.time()
	for i in range(len(dat_multi_all[1])):
		c+=1
		event_id = dat_multi_all[1][i]
		m = dat_multi_all[2][i][0].getMass()
		nstars = []
		nstars.append(len(dat_multi_all[2][i])-1)

		if error == 3:
			ev1 = [event_id, [dat_multi_all[0][i][x][0] for x in range(len(dat_multi_all[0][i]))]]
		else:
			ev1 = [event_id, dat_multi_all[0][i][0]]
		for j in range(len(dat_multi_fps[1])):
			if dat_multi_fps[1][j] == event_id:
				if error == 3: 
					ev1.append([dat_multi_fps[0][j][x][0] for x in range(len(dat_multi_fps[0][j]))])
					nstars.append(len(dat_multi_fps[2][j])-1)
				else:
					ev1.append(dat_multi_fps[0][j][0])
					nstars.append(len(dat_multi_fps[2][j])-1)

		if len(ev1) == 2: 
			ev1.append(None)
			nstars.append(0)

		for j in range(len(dat_multi_sig[1])):
			if dat_multi_sig[1][j] == event_id:
				if error == 3: 
					ev1.append([dat_multi_sig[0][j][x][0] for x in range(len(dat_multi_sig[0][j]))])
					nstars.append(len(dat_multi_sig[2][j])-1)
				else: 
					ev1.append(dat_multi_sig[0][j][0])
					nstars.append(len(dat_multi_sig[2][j])-1)

		if len(ev1) == 3: 
			ev1.append(None)
			nstars.append(0)
		acc = -999
		ind = -1
		counter =0 
		for j in range(len(dat_single[1])):
			if dat_single[1][j] == event_id:
				counter+=1
				if ind == -1:
					if error == 3: 
						acc = percentile(np.array([dat_single[0][j][x][0][6]/dat_single[0][j][x][0][2] \
								for x in range(len(dat_single[0][j]))]))[1]
					else:
						acc = dat_single[0][j][0][3]
					ind = j
				else:
					if error == 3:
						pc = percentile(np.array([dat_single[0][j][x][0][6]/dat_single[0][j][x][0][2] \
								for x in range(len(dat_single[0][j]))]))[1]
						if pc > acc:
							acc = pc
							ind = j
					else:
						if dat_single[0][j][0][3] > acc:
							acc = dat_single[0][j][0][3]
							ind = j
		if ind != -1: 
			nstars.append(1)
			if error == 3: 
				ev1.append([dat_single[0][ind][x][0] for x in range(len(dat_single[0][ind]))])
				#if percentile(np.array([dat_multi_all[0][i][x][0][6]/dat_multi_all[0][i][x][0][3] \
				#					for x in range(len(dat_multi_all[0][i]))]))[1] > -1:
			else: 
				ev1.append(dat_single[0][ind][0])
				#if dat_multi_all[0][i][0][3] > 0:
		else:
			ev1.append(None)
			nstars.append(0)
		ev1.append(nstars)


		ev1.append(m)
		ev.append(ev1)

	t[0] = time.time()-ti

	t0 = time.time()
	#select_events 
	#plot line
	if len(Folder) != 0: 
		if Folder[-1] != '/':Folder + '/'
	if not os.path.isdir(imagepath +Folder): os.mkdir(imagepath+ Folder)
	outstrings = []
	if k is None:
		k = np.arange(len(ev))
	for kk in range(len(k)):
		ti = time.time()

		fig = plt.figure(figsize = (12,11))
		ax = fig.add_subplot(111)
		plt.subplots_adjust(
				top=0.9,
				bottom=0.1,
				left=0.16,
				right=0.9,
				hspace=0.2,
				wspace=0.2)	

		i = abs(k[kk])
		ax.tick_params(axis = 'both', which = 'major', labelsize = 25)

		ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
		#if kk%2 == 0:
		ax.set_ylabel('$\Delta M \,/\, M_{int} \,[\%]$', fontsize = 30)
		ax.set_xlabel('Number of used background stars', fontsize = 30)

		ax.set_title(name(ev[i][0]), fontsize = 40)
		plt.tick_params(
			axis='x',			# changes apply to the x-axis
			which='both',		# both major and minor ticks are affected
			bottom=False,		# ticks along the bottom edge are off
			top=False,			# ticks along the top edge are off
			labelbottom=True)
		xx =[]
		num = []
		if error == 3:
				#for k in range(50): 
				ll = 0.2 

				color = iter(cc)


				t[4] = time.time()-ti	
				print(name(ev[i][0]), ev[i][-1]) 

				for j in range(1,5): 
					ti = time.time()
					sep = 0.5*j

					c = next(color)
					good = True 		
					if j <3: 
						if ev[i][-2][j-1] == ev[i][-2][j]: good = False
						
					if good:
						if ev[i][j] is not None:
								mm = np.array([ev[i][j][x][6]/ev[i][j][x][2]*100 for x in range(len(ev[i][j])) if ev[i][j][x][4] != -999])
								mp = np.array([ev[i][j][x][7]/ev[i][j][x][2]*100 for x in range(len(ev[i][j])) if ev[i][j][x][4] != -999])
								mm1 = mm[np.where(mm > -200)]
								mp1 = mp[np.where(mp < 200)]
								mm1 = mm1[np.where(mm1 < 0)]
								mp1 = mp1[np.where(mp1 > 0)]		

								parts = ax.violinplot([mp1], positions = [sep] , showextrema = False)
								for pc in parts['bodies']:
									pc.set_facecolor(c)
									pc.set_edgecolor(c)
									pc.set_alpha(0.7)
									pc.set_zorder(-0.1*j)
								rm = percentile(mm)
								rp = percentile(mp)			

								plt.plot([sep,sep, sep],rp[0:3] ,'_-',c = black, markersize = 30,linewidth = 8,markeredgewidth = 8)	
								print(j, np.array(rp[0:3]) *  ev[i][-1]/100, rp[0:3])		
								if j == 4: 
									xlim = ax.get_xlim()
									plt.plot([0,sep*5/4], [rp[1],rp[1]], linestyle= '--',dashes=(10, 5), linewidth = 5,color = black,zorder=-20)
								ylim = ax.get_xlim()
								xx.append(sep)
								num.append(str(ev[i][-2][j-1]))
								#plt.text(sep, ylim[0],str(ev[i][-1][j-1]), verticalalignment = 'bottom',\
								#		horizontalalignment = 'right', fontsize = 40)
								t[2] += time.time()-ti
						else: 
							xx.append(sep)
							num.append(' ')
					else: 
						xx.append(sep)
						num.append(' ')



		else:
			alpha = 1
			ll = 5
			color=iter(([0, 0., 1., alpha ],\
					[0, 1, 0., alpha ], \
					[1, 0, 0.5, alpha ], \
					[1., 0.5,0, alpha ]))
			for j in range(1,len(ev[i])): 
				c = next(color)
				if ev[i][j][0] is not None:
					m_minus,m,mplus = ev[i][j][3:6]
					m = ev[i][j][2]
					plt.plot([sep,sep,sep],[m_minus/m,1,mplus/m],'.', c = c, linewidth = ll)
		ax.grid(axis = 'both',linewidth = 3, zorder = -50000)
		plt.xlim([0,sep/4*5])
		ti = time.time()

		plt.xticks(xx, num)
		plt.xticks( fontsize = 25)
		plt.yticks( fontsize = 25)
		nn = ''.join(('_'.join(name(ev[i][0]).split(' '))).split(':'))
		fig.savefig(imagepath + Folder+string +'_'+nn+'.png', format = 'png')
		outstrings.append(imagepath + Folder+string +'_'+nn+'.png')
		if paperpath is not None: fig.savefig(paperpath + string +'_'+nn +'.png', format = 'png')
		print('Create Image: '+ imagepath+ Folder + string + '_'+nn  + '.png')

		plt.close(fig)
		t[5] += time.time()-ti

	t[1] = time.time()-t0

	return outstrings
	#setup axis 

def PlotMassAccuracy(stats_vary,ty = 'mass', string = '',  percent = False, ):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	creats an mass accuracy plot for good events with vary masses
 	'''

	'''
	stats_vary: (result from Analysis(vary =True)
		0D == Data, Names
		1D == events
		2D == different inputs 
		3D == different parameters 
 		4D == num, type, init, 16p, 50p,84p, sig_m,sig_p 
	'''
	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
	else: black = color_own([0.,0.,0.,1])
	if ty == 'mass' or ty =='m':
		if string =='':
			string = 'mass_accuracy_plot_mass'
		k = 0 
	if ty == 'ra':
		if string =='':
			string = 'mass_accuracy_plot_ra'
		k = 1
	if ty == 'dec':
		if string =='':
			string = 'mass_accuracy_plot_dec'
		k = 2
	if ty == 'pmra':		
		if string =='':
			string = 'mass_accuracy_plot_pmra'
		k = 3 
	if ty == 'pmdec':
		if string =='':
			string = 'mass_accuracy_plot_pmdec'
		k = 4
	if ty == 'px':
		if string =='':
			string = 'mass_accuracy_plot_px'
		k = 5
	if ty == 'c':
		if string =='':
			string = 'mass_accuracy_plot_closest_aproach'
		k = 6

	if len(stats_vary) > 1:
		stats = stats_vary[0]
		Eventnames = stats_vary[1]
		Events = stats_vary[2]
	else: 
		stats = stats_vary[0]
		Eventnames = None
	#pic = np.arange(len(Eventnames))
	#pic = np.random.choice(pic,30,replace = False)
	#print(pic)
	pic = np.array([63, 6, 46, 23, 50, 54, 39, 2, 115, 104, 4, 24, 87, 91])
	#pic =  np.array([ 128, 38, 62, 86,  24, 72, 93, 70, 37, 84, 148, 149, 93,53, 147, 31, 90, 24, 18])
	pic =  np.unique(np.concatenate((pic,pic+10), axis = 0))

	pic = pic[np.where(pic < len(Events))]
	q = np.array([stats[i][0][0][7] for i in pic])
	q2 = np.array([stats[i][0][0][2] for i in pic])
	#pic = pic[np.where(q > 0.05)]
	fig = plt.figure(figsize = (12,8))	
	ax = fig.add_subplot(111)	
	color=iter(color_own(len(pic)+2))
	
	for i in pic:
			if k != 0: 
				if k ==6:
					Obsdate, Scandir, n_ccd = RealData.getRealObs(Events[i][0], extended = True) 
					dra = np.array([j[1][2]-j[6][2] for j in stats[i]]).reshape([1,-1])
					ddec = np.array([j[2][2]-j[7][2] for j in stats[i]]).reshape([1,-1])
					dpmra = np.array([j[3][2]-j[8][2] for j in stats[i]]).reshape([1,-1])
					dpmdec = np.array([j[4][2]-j[9][2] for j in stats[i]]).reshape([1,-1])
					dpx = np.array([j[5][2]-j[10][2] for j in stats[i]]).reshape([1,-1]) 
					cd = np.array([cosdeg(j[2][2]) for j in stats[i]]).reshape([1,-1])
					tt = np.array(Obsdate).reshape([-1,1])
					#np.abs((dra**3600000*cd+tt*dpmra)*s_scandir +  (dec).reshape(-1,1) * sc_scandir[::-1])
					tmin = int(np.mean(np.argmin(np.sqrt(np.square(dra*3600000*cd+tt*dpmra)+np.square(ddec*3600000+tt*dpmdec)), axis = 0)))
					tt = tt[max([0,tmin-10]):min([len(tt)-1,tmin+10])].reshape([-1,1])
					Scandir =Scandir[max([0,tmin-10]):min([len(Scandir)-1,tmin+10])].reshape([-1,1]) 
					s_scandir = np.sin(np.array(Scandir).reshape([-1,1]))
					c_scandir = np.sin(np.array(Scandir).reshape([-1,1]))
					mass = np.array([j[0][2] for j in stats[i]]).reshape([1,-1])

					TE = Events[i][0].getMass() * dpx * const_Einsteinradius*const_Einsteinradius
					in_par = np.min(np.sqrt(np.square(dra*3600000*cd+tt*dpmra)+np.square((ddec*3600000+tt*dpmdec)))/np.sqrt(TE), axis = 0)
					#in_par = np.min(np.sqrt(np.square(dra*3600000*cd+tt*dpmra)+np.square((ddec*3600000+tt*dpmdec))), axis = 0)
					#in_par = in_par-np.min(in_par)
					#in_par = in_par/np.mean(in_par)
					if np.max(in_par)> 2000: continue
				else:
					in_par = np.array([j[k][2]-j[k+5][2]- stats[i][0][k][2]+stats[i][0][k+5][2] for j in stats[i]])



			else: in_par = np.array([j[k][2] for j in stats[i]])
			mass = np.array([j[0][2] for j in stats[i]])
			sig_m = np.array([j[0][6] for j in stats[i]])
			sig_p = np.array([j[0][7] for j in stats[i]])
			order = np.argsort(in_par)
			sig_m = sig_m[order]
			sig_p = sig_p[order]
			in_par = in_par[order]
			#if k != 0: mass = mass/mass[-1]
			if min(in_par)> 0 or k!= 0:
				c=next(color)
				a = np.where(sig_p > np.mean(sig_p))
				b = np.where(sig_p <= np.mean(sig_p))
				ta = np.argsort(in_par[a])
				tb = np.argsort(-in_par[b])
				a2 = np.where(sig_p > np.mean(sig_m))
				b2 = np.where(sig_p <= np.mean(sig_m))
				ta2 = np.argsort(in_par[a2])
				tb2 = np.argsort(-in_par[b2])
				#plt.text(mass[0],sig_p[0] , str(i), fontsize = 20,verticalalignment = 'center',horizontalalignment = 'center')
				p1 = Polygon(np.vstack([[np.hstack([in_par[a][ta],in_par[b][tb]])],[np.hstack([sig_p[a][ta],sig_p[b][tb]])]]).T, True)
				p = PatchCollection([p1,], color = c, alpha=0.7, linewidth = 4)
				ax.add_collection(p)
				if k==0 :
					if Eventnames is None:
						if percent: plt.loglog(in_par,sig_p/mass,'.', c = c, markersize = 15)
						else:plt.loglog(in_par,sig_p, '.', c = c, markersize = 15)
					else:
						if percent: plt.loglog(in_par,sig_p/mass, '.', c = c, label = Eventnames[i], markersize = 15)
						else: plt.loglog(in_par,sig_p, '.', c = c, label = Eventnames[i], markersize = 15)
				else:
					if Eventnames is None:
						if percent: plt.plot(in_par,sig_p/mass,'.', c = c, markersize = 15)
						else:plt.plot(in_par,sig_p, '.', c = c, markersize = 15)
					else:
						if percent: plt.plot(in_par,sig_p/mass, '.', c = c, label = Eventnames[i], markersize = 15)
						else: plt.plot(in_par,sig_p, '.', c = c, label = Eventnames[i], markersize = 15)

	ax.grid(axis = 'x',linewidth = 3, zorder = -50)

	xlim = np.array(plt.xlim())
	xlim = np.array([xlim[0],2.0])

	ylim = plt.ylim()
	if k ==0:
		for i in [0.15,0.3,0.5,1]:
			plt.plot(xlim,xlim*i, color = plt.rcParams['grid.color'], linewidth = 3,  zorder = -50)
			if i != 1 :
				plt.text(2.3,2.2*i , str(int(i*100))+ '%', fontsize = 25,verticalalignment = 'center',horizontalalignment = 'center')
			else: plt.text(ylim[1]*1.15,ylim[1]*1.15 , str(int(i*100))+ '%', fontsize = 25, verticalalignment = 'center',horizontalalignment = 'center')
		plt.xlim(xlim)
		plt.ylim(ylim)

		ax.tick_params(axis = 'both', which = 'major', labelsize = 25)
		ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
		plt.xticks([0.1,0.2,0.4,1,2], [0.1,0.2,0.4,1.0,2.0], fontsize = 25)
		plt.xticks(fontsize = 25)
		plt.yticks([0.02,0.05,0.1,0.2,0.4,1], [0.02,0.05,0.1,0.2,0.4,1.0], fontsize = 25)
		plt.yticks(fontsize = 25)
	if k == 6:
		plt.xlabel(r'$\delta\theta_{min}\,\, [\theta_{E}]$', fontsize = 30)
		#plt.xlim([50,100])
		#plt.ylim([0,.5])
		ax.grid(axis = 'y',linewidth = 3, zorder = -50)

	elif k != 0:
		plt.xlabel(stats[0][0][k][1], fontsize = 30)
	else: plt.xlabel('$M\,[M_{\odot}]$', fontsize = 30)
	if percent: 
		plt.ylabel('$\sigma_{M}\,[\%]$', fontsize = 30)
		fig.savefig(imagepath + string + '_percent_' + '.png')
		if paperpath is not None: fig.savefig(paperpath + string + '_percent_' +  '.png')
		print('Create Image: '+ imagepath+ string + '_percent_' + '.png')
	else: 
		plt.ylabel('$\sigma_{M}\,[M_{\odot}]$', fontsize = 30)
		fig.savefig(imagepath + string + '.png')
		if paperpath is not None: fig.savefig(paperpath + string + '.png')
		print('Create Image: '+ imagepath+ string + '.png')
	plt.close(fig)

def PlotMassInOut(stats_vary,string = 'mass_input_output'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	if string == '':
		string = 'mass_input_output'

	dark= 'black' in plt.rcParams['savefig.facecolor']
	if dark: 
		string = 'dark_'+string
		black = color_own([0.7,0.7,0.7,1])
	else: black = color_own([0.,0.,0.,1])
	gray = color_own([0.5,0.5,0.5,1])
	if isinstance(stats[0][0][0],int):
		m_minus = np.array([i[0][3] for i in stats])
		stats = [stats[j] for j in np.where(m_minus>0)[0]]
		print(len(stats))

	fig = plt.figure(figsize = (12,8))

	ax1 = plt.subplot2grid((3, 1), (0, 0),rowspan=2)
	color=iter(color_own(len(stats)+1))
	for i in range(len(stats)):
		if isinstance(stats[i][0][0],int): 
			mass = stats[i][0][2]
			out_mass = stats[i][0][4]
			sig_p_mass = stats[i][0][7]
			sig_m_mass =  -stats[i][0][6]
			if out_mass > mass: sig_mass = sig_m_mass
			else:sig_mass = sig_p_mass
		else:
			mass = np.array([j[0][2] for j in stats[i]])
			out_mass = np.array([j[0][4] for j in stats[i]])
			sig_p_mass = np.array([j[0][7] for j in stats[i]])
			sig_m_mass = np.array([-j[0][6] for j in stats[i]])
			sig_mass = sig_p_mass
			sig_mass[np.where((out_mass - mass) > 0)] = sig_m_mass[np.where((out_mass - mass) > 0)]
			order = np.argsort(mass)
			out_mass = out_mass[order]
			mass = mass[order]
			sig_mass = sig_mass[order]
		#if k != 0: mass = mass/mass[-1]
		c=next(color)
		plt.errorbar(mass,out_mass,yerr = sig_mass/np.sqrt(499) , fmt = 'o', c = 'k') 
	t = ax1.get_xlim()
	plt.plot(t,t, '-', color = gray, linewidth = 4, zorder = -50)
	plt.yticks( fontsize = 25)
	ax2 = plt.subplot2grid((3, 1), (2, 0))
	fig.subplots_adjust(hspace=0.1)  
	color=iter(color_own(len(stats)+1))
	q = []
	for i in range(len(stats)):
		if isinstance(stats[i][0][0],int):			
			mass = stats[i][0][2]
			out_mass = stats[i][0][4]
			sig_p_mass = stats[i][0][7]
			sig_m_mass =  -stats[i][0][6]
			if out_mass > mass: sig_mass = sig_p_mass
			else:sig_mass = sig_p_mass
		else:
			mass = np.array([j[0][2] for j in stats[i]])
			out_mass = np.array([j[0][4] for j in stats[i]])
			sig_p_mass = np.array([j[0][7] for j in stats[i]])
			sig_m_mass = np.array([-j[0][6] for j in stats[i]])
			sig_mass = sig_p_mass
			sig_mass[np.where((out_mass - mass) > 0)] = sig_m_mass[np.where((out_mass - mass) > 0)]
			order = np.argsort(mass)
			out_mass = out_mass[order]
			mass = mass[order]
		c=next(color)
		q.append((out_mass - mass) / sig_mass * np.sqrt(499))
		plt.plot(mass,(out_mass - mass) / sig_mass * np.sqrt(499), '.' , c = 'k') 		
	plt.plot(t,[0,0], '-', color = gray, linewidth = 4, zorder = -50)
	plt.plot(t,[-1, -1], '--', color = gray, linewidth = 3,zorder = -50)
	plt.plot(t,[1 ,1],'--', color = gray, linewidth = 3,zorder = -50)
	q = np.hstack(q)
	print(np.mean(q),np.median(q))
	ax1.set_xticklabels([])
	ax1.tick_params(axis = 'y', which = 'major', labelsize = 25,top=False, bottom=False,)
	ax1.tick_params(axis = 'y', which = 'minor', labelsize = 15,top=False, bottom=False,)
	ax2.tick_params(axis = 'both', which = 'major', labelsize = 25)
	ax2.tick_params(axis = 'both', which = 'minor', labelsize = 15)
	ax2.set_xlabel('$M_{in}\,[M_{\odot}]$', fontsize = 30)
	ax1.set_ylabel('$M_{out}\,[M_{\odot}]$', fontsize = 30)
	ax2.set_ylabel('$\Delta M/\sigma\,M$', fontsize = 30)
	t = ax1.get_xlim()


	plt.ylim([-5,5])

	fig.savefig(imagepath + string+ '.png', format = 'png')
	print('Create Image: '+ imagepath + string +'.png')
	plt.close(fig)
	
def PlotHistAll(Analysis_result, string = 'mass_accuracy_histogram', log = False, logbin = True, \
	ylim = [0,60], xlim = [0.0,1.92], extended = True, future = False, use_error = 0, vary = None, ty = 'Mass'):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	'''
	plot an histogram of input masses with color coded accuracy
	'''
	dark= 'black' in plt.rcParams['savefig.facecolor']

	if dark: string = 'dark_'+string
	try: 
		a = Analysis_result[1][0][1][1]

		TCA = None
	except IndexError:
		stats = np.array([i[0][2:] for i in Analysis_result[0]])
		if ty == 'Mass' or ty =='m':
			stats_plot = stats[:,0]
		elif ty == 'GL':
			stats_plot = np.array([i[0].getMag() for i in Analysis_result[2]])
		elif ty == 'GS':
			stats_plot = np.array([i[1].getMag() for i in Analysis_result[2]])

		elif ty == 'dG':
			stats_plot = np.array([i[1].getMag() - i[0].getMag() for i in Analysis_result[2]])
		elif ty == 'c': 
			stats_plot = np.array([dist(i[0],i[1],unit = 'arcsec', T = i[1].getTca()) for i in Analysis_result[2]])

		stats_in = stats
		stats_v = []
		stats_in_v = []
		stats_plot_v = []
		vary_bool =  vary is not None
		TCA_v = []
		if vary_bool:
			Events_v= vary[2]
			if ty == 'Mass' or ty =='m':
				stats_plot_v = np.array([i[0].getMass() for i in vary[2]])
			elif ty == 'GL':
				stats_plot_v = np.array([i[0].getMag() for i in vary[2]])
			elif ty == 'GS':
				stats_plot_v = np.array([i[1].getMag() for i in vary[2]])
			elif ty == 'dG':
				stats_plot_v = np.array([i[1].getMag() - i[0].getMag() for i in vary[2]])
			elif ty == 'c': 
				stats_plot_v = np.array([dist(i[0],i[1], unit = 'arcsec',T = i[1].getTca()) for i in vary[2]])
			for i in range(len(vary[0])):
					stats_in_v.append(Events_v[i][0].getMass())
					mm_v = np.array([vary[0][i][x][0][6] for x in range(len(vary[0][i])) if vary[0][i][x][0][4] != -999])
					mp_v = np.array([vary[0][i][x][0][7] for x in range(len(vary[0][i])) if vary[0][i][x][0][4] != -999])
					mmm_v = percentile(mm_v)
					ppp_v = percentile(mp_v)
					delta_M_v = max(ppp_v[1], -mmm_v[1])/Events_v[i][0].getMass()
					stats_v.append(delta_M_v)
					TCA_v.append(Events_v[i][0].getTca())

			stats_v = np.array(stats_v)
			stats_in_v = np.array(stats_in_v)
			stats_plot_v  = np.array(stats_plot_v)
			TCA_v = np.array(TCA_v)
			
			gmag = np.array([i[0].getMag() for i in Events_v]) 
			Good_ = np.array([name_good(i[0].getId()) for i in Events_v]) 
			gg = np.where((gmag > 5) & (Good_))
			
			stats_v = stats_v[gg]
			stats_plot_v = stats_plot_v[gg]

			TCA_v = TCA_v[gg]

			if not extended:
				ngm_v = np.where(TCA_v < 2019.5)
			if future:
				ngm_v = np.where(TCA_v > 2019.5)
				stats_v = stats_v[ngm_v]
				stats_plot_v = stats_plot_v[ngm_v]



		Eventnames =Analysis_result[1] 
		Events = Analysis_result[2] 
		ID = np.array([i[0].getId() for i in Events]) 
		gmag = np.array([i[0].getMag() for i in Events]) 
		Good_ = np.array([name_good(i[0].getId()) for i in Events]) 
		gg = np.where((gmag > 5) & (Good_))
		stats = stats[gg]
		stats_plot = stats_plot[gg]
		stats_in = stats_in[gg]
		ID = ID[gg]
		TCA = np.array([ev[0].getTca() for ev in Events])
		TCA = TCA[gg]
		if not extended:
			ngm = np.where(TCA < 2019.5)[0]
		if future:
			ngm = np.where(TCA > 2019.5)[0]
			stats = stats[ngm]
			stats_plot = stats_plot[ngm]



	else: 
		stats = np.array([i[0][2:] for i in Analysis_result])
		stats_plot = np.array([i[0][2:] for i in Analysis_result])

	
	#--------------------------------------------
	#limits for the colorcoding
	lim1 = 0.15
	lim2 = 0.000
	lim3 = 0.3
	lim4 = 0.5
	lim = [lim1, lim2, lim3, lim4, 1]

	if vary_bool:
		p10 = np.where(stats_v < lim1)[0]
		p20 = np.where(stats_v < lim2)[0]
		p30 = np.where(stats_v < lim3)[0]
		p50 = np.where(stats_v < lim4)[0]
		p100 = np.where(stats_v < 1)[0]
	elif use_error==-1:
		p10 = np.where(-stats[:,4] / stats[:,0] < lim1)[0]
		p20 = np.where(-stats[:,4] / stats[:,0] < lim2)[0]
		p30 = np.where(-stats[:,4] / stats[:,0] < lim3)[0]
		p50 = np.where(-stats[:,4] / stats[:,0] < lim4)[0]
		p100 = np.where(-stats[:,4] / stats[:,0] < 1)[0]
	elif use_error==1:
		p10 = np.where(stats[:,5] / stats[:,0] < lim1)[0]
		p20 = np.where(stats[:,5] / stats[:,0] < lim2)[0]
		p30 = np.where(stats[:,5] / stats[:,0] < lim3)[0]
		p50 = np.where(stats[:,5] / stats[:,0] < lim4)[0]
		p100 = np.where(stats[:,5] / stats[:,0] < 1)[0]
	elif use_error == -999:
		p10 = np.where(np.minimum(stats[:,5],-stats[:,4]) / stats[:,0] < lim1)[0]
		p20 = np.where(np.minimum(stats[:,5],-stats[:,4]) / stats[:,0] < lim2)[0]
		p30 = np.where(np.minimum(stats[:,5],-stats[:,4]) / stats[:,0] < lim3)[0]
		p50 = np.where(np.minimum(stats[:,5],-stats[:,4]) / stats[:,0] < lim4)[0]
		p100 = np.where(np.minimum(stats[:,5],-stats[:,4]) / stats[:,0] < 1)[0]
	else:
		p10 = np.where(np.maximum(stats[:,5],-stats[:,4]) / stats[:,0] < lim1)[0]
		p20 = np.where(np.maximum(stats[:,5],-stats[:,4]) / stats[:,0] < lim2)[0]
		p30 = np.where(np.maximum(stats[:,5],-stats[:,4]) / stats[:,0] < lim3)[0]
		p50 = np.where(np.maximum(stats[:,5],-stats[:,4]) / stats[:,0] < lim4)[0]
		p100 = np.where(np.maximum(stats[:,5],-stats[:,4]) / stats[:,0] < 1)[0]
	p = [p10,p20, p30, p50, p100]

	#--------------------------------------------
	
	
	if not extended and TCA is not None:
		print(len(p10),len(p20),len(p30),len(p50),len(p100), len(stats), len(ngm))
	else:
		print(len(p10),len(p20),len(p30),len(p50),len(p100), len(stats), len(np.unique(ID)))

	#--------------------------------------------
	#setup figure
	fig = plt.figure(figsize = (12,8))	
	ax = fig.add_subplot(111)
	plt.subplots_adjust(
			top=0.9,
			bottom=0.16)	
	#--------------------------------------------
	
	#--------------------------------------------
	#defines a scale for the bins
	if ty == 'Mass' or ty == 'm' : 
		if log:
			plt.xscale('log')
			bins = 1.3 ** np.arange(0,20) * 0.05 
		elif logbin: 
			d = 1.25 ** np.arange(0,11) * 0.05
			bins = [np.sum(d[:i])+0.05 for i in range(len(d))]
		else:
			bins = 0.05 + 0.05 * np.arange(0,30)
	elif ty == 'c':
		bins = [0,1/3,2/3,1,4/3,5/3,2,3,4,5,6,7]
		#bins = 0.5*np.arange(14)
		#bins = 0.35 * 20 ** (np.arange(-2,5)/4)
		print(bins)
	elif 'G' in ty:
		bins = 8 + np.arange(0,11)
	else:
		bins = 10
	#--------------------------------------------
	
	#--------------------------------------------
	#setup colores and hatches
	col_list = [[0.6,0.6,0.6,0.6],[1,0,0,1],[1,1,0,1],[0,0,1,1],[0,1,0,1]]
	color=iter(color_own(col_list))

	hatch = iter(['', '//', 'x', '\\', '|'])
	if dark:
		black = color_own([0.7,0.7,0.7,1])
	else: black = color_own([0.,0.,0.,1])
	#--------------------------------------------
	
	#--------------------------------------------
	# plot input distribution
	if ty == 'Mass' or ty == 'm':
		if not extended and TCA is not None:	
			plt.hist(stats_in[ngm,0], bins = bins, ls = '-', histtype = 'step',linewidth = 2, color = black, label = '5 years sample', zorder = -2)
		if future and TCA is not None:
			q = plt.hist(stats_in[ngm,0], bins = bins, histtype = 'step', color = black,linewidth = 5, label = 'future events', zorder = 55)
		else:
			q = plt.hist(stats_in[:,0], bins = bins, histtype = 'step', color = black,linewidth = 5, label = '10 years sample', zorder = 55)
	else:
		q = plt.hist(stats_plot, bins = bins, histtype = 'step', color = black,linewidth = 5, label = '10 years sample', zorder = 55)

	#--------------------------------------------


	#--------------------------------------------
	# plot distribution for each limit
	for j in range(4, -1, -1): 
		c = next(color)
		h = next(hatch)
		if len(p[j]>0):
			if vary_bool:
				print(len(p[j]))
				plt.hist(stats_plot_v[p[j]], bins = bins, color = c, hatch = h, rwidth = 0.9, \
					label = '$\sigma < %i\%s$'%(100*lim[j],'%'))
			else:
				plt.hist(stats_plot[p[j]], bins = bins, color = c, hatch = h, rwidth = 0.9, \
					label = '$\sigma < %i\%s$'%(100*lim[j],'%'))
	#--------------------------------------------	

	#--------------------------------------------
	#scale y axis with labels 
	if max(q[0]) <  ylim[1]:
		ylim = plt.ylim()
	n = np.where(q[0] > ylim[1])[0]

	if ylim[1] > 20:
		y_minor_ticks = np.arange(ylim[0], ylim[1]+5,1)
		y_major_ticks = np.arange(ylim[0]+10, ylim[1]+10,10)
	else:
		y_minor_ticks = np.arange(ylim[0], ylim[1]+5,1)
		y_major_ticks = np.arange(ylim[0]+2, ylim[1]+5,2)
	#--------------------------------------------

	#--------------------------------------------
	#setup axis and labels 
	ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
	ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
	ax.set_yticks(y_major_ticks)
	ax.set_yticks(y_minor_ticks, minor = True)
	if 'G' in ty: plt.xticks([(round(b,2))  for b in bins],["%i" % (round(b,2))  for b in bins],  fontsize = 20)
	elif not isinstance(bins,int):plt.xticks([(round(b,2))  for b in bins],["%0.2f" % (round(b,2))  for b in bins],  fontsize = 20)
	for tick in ax.get_xticklabels():
		tick.set_rotation(90)
	plt.yticks( fontsize = 25)
	

	ax.set_axisbelow(True)
	ax.yaxis.grid(color='gray', linestyle='dashed', linewidth = 3, zorder = -20)
	for j in n:
		plt.text((q[1][j]+q[1][j+1])/2, ylim[1]-2,str(int(q[0][j])), verticalalignment = 'top',\
			horizontalalignment = 'center', fontsize = 16)
	if ty == 'Mass' or ty == 'm': plt.xlim(xlim)
	plt.ylim(ylim)
	if 'c' in ty: 
		plt.xlabel('$\phi_{min}\,$["]',fontsize = 30)
		plt.legend(fontsize = 25,loc =1)
		plt.title('Impact parameter', fontsize = 40)
	elif 'G' in ty: 
		plt.xlabel('$G\,[mag]$', fontsize = 30)
		plt.legend(fontsize = 25,loc =2)
		plt.title('G magnitude (source)', fontsize = 40)
	elif ty == 'Mass'or ty =='m': 
		plt.xlabel('$M\,[M_{\odot}]$', fontsize = 30)
		plt.legend(fontsize = 25,loc =1)
		if future: plt.title('Future events', fontsize = 40)
		elif extended: plt.title('Extended mission', fontsize = 40)
		else : plt.title('Nominal mission', fontsize = 40)
	else:plt.legend(fontsize = 25,loc =0)
	plt.ylabel('#', fontsize = 30)

	
	
	#--------------------------------------------

	#--------------------------------------------
	#save figure
	if ty != 'Mass' and ty != 'm':
		string = string+ '_'+ty
	fig.savefig(imagepath + string +'.png', format = 'png')
	if paperpath is not None: fig.savefig(paperpath + string + '.png', format = 'png')
	print('Create Image: '+ imagepath+ string + '.png')
	plt.close(fig)
	#--------------------------------------------

def PlotCHI2(MC_res10):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	nn = []
	cc = 0
	for i in range(len(MC_res10['Results'])):
		fig = plt.figure(figsize = (12,10))
		plt.ylabel('#', fontsize = 26)
		plt.xlabel('$\chi^2_{reduced}$', fontsize = 26)
		plt.xticks( fontsize = 20)
		plt.yticks( fontsize = 20)
		data = MC_res10['Chi2'][i]
		if isinstance(data,int): pass
		else:
			plt.hist(data[np.where(data<10)])
			name = str(MC_res10['Eventlist'][i][0].getId())
			if name in nn: 
				ii = 1
				while name+'_'+str(ii) in nn:
					ii+=1
				name = name+'_'+str(ii) 
			nn.append(name)
			plt.title(name,  fontsize = 30)
			cc+=1
			fig.savefig(imagepath + 'CHI2/CHI2_DIST_' + name + '.png')
			plt.close(fig)
	print('Create %i Images: '%(cc)+ 'CHI2_DIST')
