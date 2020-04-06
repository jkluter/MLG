__all__ = ['sindeg','cosdeg','sinmas','cosmas','dist','dist_scandir','unitFac','percentile']
import  numpy as np
from MLG.const import *

def sindeg(x):
	#sin(x) for x in degree
	return np.sin(x*DEG)
def cosdeg(x):
	#cos(x) for x in degree
	return np.cos(x*DEG)
def sinmas(x):
	#sin(x) for x in mas
	return np.sin(x*MAS)
def cosmas(x):
	#cos(x) for x in mas
	return np.cos(x*MAS)


def dist(ra,dec,ra2=None,dec2=None, unit = 'deg', T = None):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	if isinstance(ra, (int, float, np.ndarray)): 
		unitfactor=1
		if 'mas' in unit: unitfactor = 3.6e6
		if 'arcsec' in unit: unitfactor = 3.6e3

		if ra2 is None:
			return ((ra[0]-dec[0])**2*cosdeg(ra[1])**2 + (ra[1]-dec[1])**2)**0.5*unitfactor
		else:
			return ((ra-ra2)**2*cosdeg(dec)**2 + (dec-dec2)**2)**0.5*unitfactor
	else: 
		try:
			star1 = ra
			star2 = dec
			radec1 = np.array(star1.getPos())
		except: 
			print(type(ra))
			print('expect Number or np.array or STAR_Objekt')
			raise TypeError

		unitfactor=1
		if 'mas' in unit: unitfactor = 3.6e6
		if 'arcsec' in unit: unitfactor = 3.6e3
		if T is None:  T = 0.
		elif isinstance(T,str): T = star1.getTca()-2015.5
		elif T > 2000: T = T-2015.5
		dpos = (np.array(star2.getPos()) - radec1)*unitfactor* np.array([cosdeg(radec1[1]),1])
		dpm_T = T * (np.array(star2.getPm(unit))-np.array(star1.getPm(unit)))
		return(np.sqrt(np.sum((dpos+dpm_T)*(dpos+dpm_T)))) 


def dist_scandir(scandir,ra,dec,ra2=None,dec2=None, unit = 'deg', ):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	# calculate the separation in along scan und cross scan direction
	sc_scandir = [np.sin(scandir),np.cos(scandir)]
	unitfactor = 1
	if 'mas' in unit: unitfactor = 3.6e6
	if 'arcsec' in unit: unitfactor = 3.6e3
	if ra2 is None:
		return np.abs((ra[0]-dec[:,0]).reshape(-1,1)*cosdeg(ra[1])*sc_scandir+ (ra[1]-dec[:,1]).reshape(-1,1) * sc_scandir[::-1]*np.array([1,-1]))*unitfactor
	else:
		return np.abs((ra-ra2).reshape(-1,1)*cosdeg(dec)*sc_scandir +  (dec-dec2).reshape(-1,1) * sc_scandir[::-1] *np.array([1,-1]))*unitfactor

def unitFac(unit_in, unit_out):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	if 'mas' in unit_in: cin = 1.e0
	if 'arcsec' in unit_in: cin	= 1.e3
	if 'deg' in unit_in: cin = 3.6e6
	if 'mas' in unit_out: cout = 1.e0
	if 'arcsec' in unit_out: cout	= 1.e3
	if 'deg' in unit_out: cout = 3.6e6
	return cin/cout


def percentile(vec):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	vec  = vec[np.where(vec != -999)]

	n = len(vec)
	if n  == 0 : return -999, -999, -999, -999, -999
	n16 = round(n/100*15.866)
	n50 = round(n/100*50)
	n84 = round(n/100*84.134)

	vec_ord = np.sort(vec)
	p16 = vec_ord[n16]
	p50 = vec_ord[n50]
	p84 = vec_ord[n84]

	sminus = (p84 - p50)
	splus = (p16 - p50)
	return p16, p50, p84, splus,sminus

