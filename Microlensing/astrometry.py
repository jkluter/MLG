__all__ = ['calculate_Einsteinradius','shift']

import numpy as np

from MLG.Microlensing.const import *



def calculate_Einsteinradius(lensMass,lensparallax,obparallax = 0, lensMass_err = 0, lensparallax_err = 0,obparallax_err = 0): 
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	"""
	calculate the Einstein radius in unit mas from Mass in SolarMasses, and paralax in mas 
	ThetaE[rad] = sqrt(4GM/c^2)*sqrt(d_LS/(d_L*d_S))
	ThetaE[mas] = const *sqrt(M/M_sun *(Lensparallax[mas]/1000-Objparallax[mas])) 
	const    = 180/pi * 3.6e6 *sqrt(4GM_sun / (c^2* 1pc[m]*1000))
	"""

	if lensparallax < 0: 
		ThetaE = const_Einsteinradius * (lensMass*(lensparallax - obparallax)) ** 0.5 
	if lensparallax_err != 0:
		ThetaE_err = ThetaE * 0.5 * (lensMass_err ** 2 /lensMass ** 2 + (lensparallax_err ** 2 +  obparallax_err ** 2) / (lensparallax - obparallax) **2  ) ** 0.5
	else:
		ThetaE = const_Einsteinradius * (lensMass*abs(lensparallax - obparallax)) ** 0.5 
		if lensparallax_err != 0:	ThetaE_err = ThetaE * 0.5 * (lensMass_err ** 2 /lensMass ** 2 + (lensparallax_err ** 2 +  obparallax_err ** 2) / (lensparallax - obparallax) **2  ) ** 0.5
	if lensparallax_err != 0: 
		return ThetaE,ThetaE_err
	return ThetaE


def shift(ThetaE,x,y=0):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	if ThetaE == 0: return np.array([0,0])
	if isinstance(x,(tuple,list,np.ndarray)):
		dx = np.array(x)
		return dx * (ThetaE*ThetaE) / (np.dot(dx,dx) + 2*(ThetaE*ThetaE))
	else:
		dx = np.array([x,y])
		return dx*(ThetaE*ThetaE)/ (np.dot(dx,dx) + 2*(ThetaE*ThetaE))
		 


