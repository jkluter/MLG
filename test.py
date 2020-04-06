import numpy as np
import MLG
import random
import matplotlib.pyplot as plt
from MLG.Math import cosdeg,sindeg,percentile
import time
import importlib
from MLG.Simulation.RealData import loadRealDataMulti,loadRealData,getRealObs
from  MLG.Simulation import RealData ,RawData
from MLG.Modeling.fitting_micro import Fit_Micro
import MLG.Analysis.utils

from joblib import Parallel, delayed
import warnings


def re(test = True): 
	for i in range(2):
		if test: 
			try: importlib.reload(MLG.test)
			except: pass
		importlib.reload(MLG)
		importlib.reload(MLG.utils)

		importlib.reload(MLG.Math.utils)
		importlib.reload(MLG.Math)
		importlib.reload(MLG.name)
		#importlib.reload(MLG.image)
		importlib.reload(MLG.Analysis)
		importlib.reload(MLG.Analysis.Table)
		importlib.reload(MLG.Analysis.Plots)
		importlib.reload(MLG.Analysis.utils)
		importlib.reload(MLG.Simulation)
		importlib.reload(MLG.Simulation.MonteCarlo)
		importlib.reload(MLG.Simulation.RawData)
		importlib.reload(MLG.Simulation.RealData)

		importlib.reload(MLG.StellarMotion)
		importlib.reload(MLG.StellarMotion.motion_stars)
		importlib.reload(MLG.Modeling.fitting_micro)
		importlib.reload(MLG.Modeling.fitting_motion)
		importlib.reload(MLG.StellarMotion)
		importlib.reload(MLG)
		importlib.reload(MLG.image)

		if test: 
			importlib.reload(MLG.test)

def test_ext(string='/MC_single_test_ext_29-01-2020_28_10_500_1.pkl'):
	mc = MLG.Simulation.MonteCarlo.loadPkl(string)
	stat_reg = MLG.Analysis.statistics(mc, res_string = 'Results_comp')
	stat_ext = MLG.Analysis.statistics(mc)
	dist_tca = []
	t = [i[0].getTca() for i in mc['Eventlist']]
	for i in mc['Eventlist']:
		radec_lens =np.array(i[0].getPos()) 
		pmradec_lens = np.array(i[0].getPm()) * np.array([1/MLG.Math.cosdeg(radec_lens[1]),1])
		radec_source =np.array(i[1].getPos()) 
		pmradec_source = np.array(i[1].getPm()) * np.array([1/MLG.Math.cosdeg(radec_lens[1]),1])
		dist = lambda tobs: MLG.Math.dist(radec_lens + tobs*pmradec_lens/3.6e6, radec_source + tobs*pmradec_source/3.6e6, unit = 'mas')
		dist_tca.append(dist(i[0].getTca()-2015.5))	


	mass = np.array([i[0][2] for i in stat_reg[0]])
	mass_reg = np.array([i[0][4] for i in stat_reg[0]])
	mass_ext = np.array([i[0][4] for i in stat_ext[0]])
	sig_reg = np.array([i[0][7] for i in stat_reg[0]])
	sig_ext = np.array([i[0][7] for i in stat_ext[0]])
	return np.array([sig_reg,sig_ext,mass_reg,mass_ext,mass, t, dist_tca]).T


def test():
	mc = MLG.Simulation.MonteCarlo.loadPkl('MC_single_external_obs_21-02-2020_511_10_500_1.pkl')
	print(mc['Header'])
	ec = MLG.Analysis.statistics(mc,ext_obs = True)
	nc = MLG.Analysis.statistics(mc,ext_obs = False)
	kkk=np.array([[nc[0][j][0][8],ec[0][j][0][8]] for j in range(len(nc[0]))])
	vvv=np.array([[nc[0][j][0][7],ec[0][j][0][7]] for j in range(len(nc[0]))])
	ev = np.array(nc[2])
	ev = ev[np.where((kkk[:,0]<1)&(kkk[:,0]>0))]
	vvv = vvv[np.where((kkk[:,0]<1)&(kkk[:,0]>0))]
	kkk = kkk[np.where((kkk[:,0]<1)&(kkk[:,0]>0))]
	gmag = [k[1].getMag() for k in ev]
	tca = [k[0].getTca() for k in ev]

	dist = [MLG.Math.dist(*k,unit = 'mas', T = '') for k in ev]
	import matplotlib.pyplot as plt
	fig = plt.figure()
	plt.plot(dist,(kkk[:,0]-kkk[:,1])/kkk[:,0],'.')
	plt.title('delta_rel vs dist')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_vs_dist.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(vvv[:,0]-vvv[:,1])/vvv[:,0],'.')
	plt.title('delta vs dist')
	fig.savefig('/Users/jonas/Desktop/img/delta_vs_dist.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(kkk[:,0]-kkk[:,1]),'.')
	plt.title('delta_rel_abs vs dist')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_abs_vs_dist.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(vvv[:,0]-vvv[:,1]),'.')
	plt.title('delta_abs vs dist')
	fig.savefig('/Users/jonas/Desktop/img/delta_abs_vs_dist.png')
	plt.close(fig)
	#-----------------------------------------------------------------
	fig = plt.figure()
	plt.plot(kkk[:,0],(kkk[:,0]-kkk[:,1])/kkk[:,0],'.')
	plt.title('delta_rel vs value')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_vs_value.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(vvv[:,0],(vvv[:,0]-vvv[:,1])/vvv[:,0],'.')
	plt.title('delta vs value')
	fig.savefig('/Users/jonas/Desktop/img/delta_vs_value.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(kkk[:,0],(kkk[:,0]-kkk[:,1]),'.')
	plt.title('delta_rel_abs vs value')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_abs_vs_value.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(vvv[:,0],(vvv[:,0]-vvv[:,1]),'.')
	plt.title('delta_abs vs value')
	fig.savefig('/Users/jonas/Desktop/img/delta_abs_vs_value.png')
	plt.close(fig)
	#-----------------------------------------------------------------
	fig = plt.figure()
	plt.plot(gmag,(kkk[:,0]-kkk[:,1])/kkk[:,0],'.')
	plt.title('delta_rel vs gmag')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_vs_gmag.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(gmag,(vvv[:,0]-vvv[:,1])/vvv[:,0],'.')
	plt.title('delta vs gmag')
	fig.savefig('/Users/jonas/Desktop/img/delta_vs_gmag.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(gmag,(kkk[:,0]-kkk[:,1]),'.')
	plt.title('delta_rel_abs vs gmag')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_abs_vs_gmag.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(gmag,(vvv[:,0]-vvv[:,1]),'.')
	plt.title('delta_abs vs gmag')
	fig.savefig('/Users/jonas/Desktop/img/delta_abs_vs_gmag.png')
	plt.close(fig)
	#-----------------------------------------------------------------
	fig = plt.figure()
	plt.plot(tca,(kkk[:,0]-kkk[:,1])/kkk[:,0],'.')
	plt.title('delta_rel vs tca')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_vs_tca.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(tca,(vvv[:,0]-vvv[:,1])/vvv[:,0],'.')
	plt.title('delta vs tca')
	fig.savefig('/Users/jonas/Desktop/img/delta_vs_tca.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(tca,(kkk[:,0]-kkk[:,1]),'.')
	plt.title('delta_rel_abs vs tca')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_abs_vs_tca.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(tca,(vvv[:,0]-vvv[:,1]),'.')
	plt.title('delta_abs vs tca')
	fig.savefig('/Users/jonas/Desktop/img/delta_abs_vs_tca.png')
	plt.close(fig)
	#-----------------------------------------------------------------
	tca = np.array(tca)
	dist= np.array(dist)
	kkk=kkk[np.where(tca> 2019.5)]
	vvv=vvv[np.where(tca> 2019.5)]
	dist = dist[np.where(tca> 2019.5)]

	fig = plt.figure()
	plt.plot(dist,(kkk[:,0]-kkk[:,1])/kkk[:,0],'.')
	plt.title('delta_rel vs dist (future)')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_vs_dist_future.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(vvv[:,0]-vvv[:,1])/vvv[:,0],'.')
	plt.title('delta vs dist (future)')
	fig.savefig('/Users/jonas/Desktop/img/delta_vs_dist_future.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(kkk[:,0]-kkk[:,1]),'.')
	plt.title('delta_rel_abs vs dist (future)')
	fig.savefig('/Users/jonas/Desktop/img/delta_rel_abs_vs_dist_future.png')
	plt.close(fig)
	
	fig = plt.figure()
	plt.plot(dist,(vvv[:,0]-vvv[:,1]),'.')
	plt.title('delta_abs vs dist (future)')
	fig.savefig('/Users/jonas/Desktop/img/delta_abs_vs_dist_future.png')
	plt.close(fig)

	fig = plt.figure()
	plt.plot(vvv[:,0],vvv[:,1],'.')
	plt.title('external vs Gaia future')
	fig.savefig('/Users/jonas/Desktop/img/value_value_future.png')
	plt.close(fig)






