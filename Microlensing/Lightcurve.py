from astropy.table import Table
from astropy.time import Time

import numpy as np
from MLG.StellarMotion import stars_position, stars_position_micro, getSun
from MLG.Math import *
from MLG.utils import color_own
from MLG.image import dark

import matplotlib.pyplot  as plt
Data = Table.read('/Users/jonas/amlensing/Final_DATA/events.3phot.vot', format = 'votable')
obs1 = Table.read('/Users/jonas/Desktop/cpt/Cand2/light_587', format = 'ascii')
obs = Table.read('/Users/jonas/Desktop/cpt/Cand2_cpt/light_628', format = 'ascii')
obs2 = Table.read('/Users/jonas/Desktop/GAIACAND2_cpt-doma-1m0-10-fa16_ip.t', format = 'ascii')

obs = obs[np.where(obs['col3']< 0.05)]
obs1 = obs1[np.where(obs1['col3']< 0.05)]

t_0 = 2015.5
t_mai = Data[1]['date']

phot = ['phot_g_mean_mag','ob_phot_g_mean_mag' ]
gmag = np.array([Data[1][i] for i in phot])
gmag2 = [12.19, 13.6]
par_col= ['mass','ra', 'dec', 'pmra', 'pmdec', 'parallax', 'ob_ra', 'ob_dec', 'ob_pmra', 'ob_pmdec', 'ob_parallax']
par = np.array([Data[1][i] for i in par_col])
print(par)
par_multi= np.vstack([par,par])
epoch = np.linspace(t_mai-t_0 - 40/365, t_mai-t_0 + 40/365, 200)
epoch_vec = np.repeat(epoch.reshape(-1,1),2,axis = 1 )
num = np.arange(2)
scsc = sindeg(par[1]), cosdeg(par[1]), sindeg(par[2]),\
					max(cosdeg(par[2]), 1e-16)

earth = getSun(t_0 + epoch).T

loc = np.repeat(np.stack([(scsc[0] * earth[:,0] - scsc[1] *  earth[:,1])\
			/ max(scsc[3],1e-6),\
			scsc[1] * scsc[2] *  earth[:,0] + scsc[0] * scsc[2]\
			*  earth[:,1] - scsc[3] *  earth[:,2]]).T.reshape(-1,1,2),2,axis = 1 )
num_vec = np.repeat(num.reshape(1,-1),len(epoch), axis = 0)  
pos_vec  = stars_position_micro(par_multi,num_vec.reshape(-1),epoch_vec.reshape(-1),loc_vec = loc.reshape(-1,2),scsc = scsc, ml = True, pml = True)
A = pos_vec[1][0]
A = A[np.where(A !=0)]
f_L = 10**(-0.4*(gmag[0]))
f_S = 10**(-0.4*(gmag[1]))
f_L2 = 10**(-0.4*(gmag2[0]))
f_S2 = 10**(-0.4*(gmag2[1]))
f_L/f_S
m1 = -2.5*np.log10(f_L+f_S*A)
m2 = -2.5*np.log10(f_L+f_S)
m12 = -2.5*np.log10(f_L2+f_S2*A)

m22 = -2.5*np.log10(f_L2+f_S2)
d = 1
if d: 
	dark(1)
	c1 =  color_own([1,0.5,0,1])
	c2 = color_own([1,1,0,1])
	c3 = color_own([0,1,1,1])
	c4 = color_own([0.6,1.2,0,1])
	c5 = color_own([1,1,1,1])
else: 
	c1 =  color_own([1,0.5,0,1])
	c2 = color_own([1,0,0,1])
	c3 = color_own([0,0,1,1])
	c4 = color_own([0,1,0,1])
	c5 = color_own([0,0,0,1])


ee = Time(epoch+t_0, format = 'jyear').jd
plt.plot(ee-2450000,m12-m22+np.median(obs['col2']), color = c1, label = 'H fluxratio')
plt.plot(ee-2450000,m1-m2+np.median(obs['col2']), color = c2,  label = 'G fluxratio')
plt.errorbar(obs2['col2']-2450000,obs2['col7']- np.median(obs2['col7'])+ np.median(obs['col2']),obs2['col8'], fmt = '.', color = c4, label = 'Obs1')

#plt.errorbar(obs1['col1']-2450000,obs1['col2']- np.median(obs1['col2'])+ np.median(obs['col2']),obs1['col3'], fmt = '.', color = c4, label = 'Obs1')
plt.errorbar(obs['col1']-2450000,obs['col2'],obs['col3'], fmt = '.', color = c3, label = 'cpt')

plt.plot([8632,8632],plt.ylim(), color = c5)
plt.ylabel('Magnitude')
plt.xlabel('HJD - 2450000')
print(np.mean(obs['col3']))
plt.gca().invert_yaxis()
plt.legend()
plt.show()

