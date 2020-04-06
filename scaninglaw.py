import numpy as np
from astropy.coordinates import get_sun as get_sun_astropy
from astropy.time import Time
from astropy.table import vstack, Table
import time
path = __file__[:-13]
def trigmul(x,imax):
	#calculate sin(k*x) and cos(k*x) for k =(0 ... imax)
	k = np.arange(imax + 1)
	sk= np.sin(k*x)
	ck = np.cos(k*x)
	return ck,sk

def trigsum(ck1, sk1, ck2, sk2):
	#calculate sin(a*x1 + b*x1) and cos(a*x1 + b*x1)
	# a = -2... imax1, = -2... imax2
	# sk1,sk2 = sin(k*x) for k = 0... 
	# ck1,ck2 = cos(k*x) for k = 0...  
	l1 = len(ck1)+2
	l2 = len(ck2)+2
	ck1 = np.repeat(np.append(ck1,ck1[2:0:-1]).reshape(-1,1),l2, axis = 1)
	sk1 = np.repeat(np.append(sk1,-sk1[2:0:-1]).reshape(-1,1),l2, axis = 1)
	ck2 = np.repeat(np.append(ck2,ck2[2:0:-1]).reshape(1,-1),l1, axis = 0)
	sk2 = np.repeat(np.append(sk2,-sk2[2:0:-1]).reshape(1,-1),l1, axis = 0)
	cc = ck1 * ck2 - sk1 * sk2
	ss = sk1 * ck2 + ck1 * sk2

	return cc,ss

def getSun_long(x):
	#determine the longitude of the sun for x in julian date
	#for vectorised computing
	return get_sun_astropy(Time(x, format = 'jd')).ra.value
getSun_long_vec = np.vectorize(getSun_long)

def sun_long(xjd, sun_ra = None):
	''' ------------------------------------------------------------------------------------------------------
	determine the longitude of the sun for x in julian date
	xjd: 	Julian date
	sun_ra:	Tabel of dates and longitudes for interpolation
	------------------------------------------------------------------------------------------------------'''
	if sun_ra != None:
		index = np.where(sun_ra[0]>xjd)[0]
		if len(index) > 0:
			index = index[0]
			if index != 0:
				return sun_ra[1][index-1] + (sun_ra[1][index] - sun_ra[1][index-1])/ (sun_ra[0][index]-sun_ra[0][index-1]) * (xjd-sun_ra[0][index-1])
	return get_sun_astropy(Time(xjd, format = 'jd')).ra.value

def dsl_dt(xjd, sun_ra = None):
	if sun_ra != None:
		index = np.where(sun_ra[0]>xjd)[0]
		if len(index) > 0:
			index = index[0]
			if index != 0:
				return (sun_ra[1][index] - sun_ra[1][index-1])/ (sun_ra[0][index]-sun_ra[0][index-1])
	dt = 0.000114155
	return  (get_sun_astropy(Time(xjd, format = 'jd')).ra.value \
			- get_sun_astropy(Time(xjd-dt, format = 'jd')).ra.value)/dt

def scaninglaw(date1=2457206.375,date2=0, xjd_ref=2457206.375 , xi = 45, xk=5.8, etaz=0, phiz=0, \
	omegaz = 59.960499070315, csba=(0.5983246005706591,0.8012538126910607),\
	parxi = None, sun_ra = None):
	
	#compunting-time-tracker
	cputt = np.zeros(10)
	''' ------------------------------------------------------------------------------------------------------
		INPUT
			date1   : date1+date2  = date in   julian days  (see note below)
			date2   :
			xjd_ref : date in julian days to which the initial conditions are referred (nominally 01 Sep 2013 TCB)
			xi      : sun aspect angle (obliquity) in degrees (nominally 45°)
			xk      : precession frequency (number of revolutions per year: denoted 'K' in the reference papers)
						- Could be < 0 for reverse SL
			etaz    : initial phase of the precession motion  in deg at xjd_ref. (\alpha in the reference papers)
			phiz    : initial phase of the rotational motion (heliotropic convention here) in deg at xjd_ref. 
					  (\Omega_0 in the reference papers)
			omegaz  : nominal scanning speed in deg/hour (= "/s  ~ 60.0)
			csba	: cos, sin of the half of the basic angle (default = (cos(53.25°), sin(53.25°)))
			parxi	: precalculated parameters for xi 
			sun_ra  : precalculated longitudes of the sun
	------------------------------------------------------------------------------------------------------
		OUTPUT  : values of parameters at date1 + date2
			xlsun   : ecliptic longitude of the sun in degrees (conventional expression of low accuracy for the Sl)
			eta     : precesion angle (v in the reference paper)
			psi     : ascending node wrt sun of the scanning circle in deg in the ecliptic 
						(not from gamma, but from Sun)
			teta    : inclination of the scanning circle in deg wrt eclitpic
			omega   : heliotropic phase of the x-axis (midway between the two fovs) in deg (\Omega in the papers)
			phi     : third euler angle about the spin axis for the x-axis in deg wrt node with ecliptic
			etap    : d(eta)/dt in "/s (= deg/h)
			psip    : dpsi/dt   in "/s (= deg/h) (including the solar motion : inertial psip  :: very bad notation !!)
			tetap   : dteta/dt  in "/s (= deg/h)
			omegap  : domega/dt in "/s (= deg/h)
			phip    : dphi/dt   in "/s (= deg/h)
			omenode : heliotropic abscissa of the node (Omega of the node. This is also  omega - phi)
			date    : epoch == date1+date2	
			z       : cartesian vector of the z axis
			pfov    : cartesian coordinates of the preceding field of view
			ffov    : cartesian coordinates of the following field of view
			pfov_dir: dirrection vector of the preceding field of view
			ffov_dir: dirrection vector of the following field of view

			cputt   : tracker of the computation time	
	------------------------------------------------------------------------------------------------------
		Note  for DATE1, DATE2 with SOFA conventions :
			The epoch DATE1 + DATE2 is a Julian Date, apportioned in any
			convenient way between the arguments date1 and date2.  For
			example, JD  =2450123.7 could be expressed in any of these
			ways, among others:
		
				  DATE1          DATE2
				
				2450123.7D0        0D0        (JD method)
				2451545D0      -1421.3D0     (J2000 method)
				2400000.5D0     50123.2D0     (MJD method)
				2450123.5D0       0.2D0       (date & time method)
		 
			The JD method is the most natural and convenient to use in
			cases where the loss of several decimal digits of resolution
			is acceptable.  The J2000 method is best matched to the way
			the argument is handled internally and will deliver the
			optimum resolution.  The MJD method and the date & time methods
			are both good compromises between resolution and convenience.
	------------------------------------------------------------------------------------------------------'''

	#----------------------------------------------------------------------------------
	#constans
	#
	raddeg = 180/np.pi
	asdeg = 1/3.6e6
	ro = 360.0
	year = 365.25
	dayh = 24.
	imax = 8

	# ----------------------------------------------------------------------------------

	#----------------------------------------------------------------------------------
	#	Calculate parameters
	#
	if parxi == None: first = True
	else:

		try:
			xiref, cxi, sxi, cxi2, sxi2, cxi3, sxi3, cxi4, sxi4, cxi5, sxi5, cxi6, sxi6, cxi8,\
				sxi8, cxi10, sxi10, cxi12, sxi12, psk2,\
				psk4, psk6, psk8, psk10, psk12, p11p1, p10p1, p21p0, p21p2, p22p2, p20p2, p31m1,\
				p31p1, p31p3, p32p1, p32p3, p33p3, p30p1, p30p3, p41p0, p41m2, p41p2, p41p4, p42p0,\
				p42p2, p42p4, p43p2, p43p4, p44p4, p40p2, p40p4,\
				qvs2, qvs4, qvs6, qvs8, q1c1, q2s2, q3c1, q3c3, q4s2, q4s4, q5c1, q5c3, q5c5,\
				q6s2, q6s4, q6s6, q7c1, q7c3, q7c5, q7c7, saa = parxi
		except: 
			first = True
		else:
			if xi == xiref: 
				first = False
			else: first = True 
	if first: 
		cputt[0] -= time.time()# ctt
		# initialisations at the first call  Disabled 10/03/2008
		# allow to activate change in the free SL parameters within a run
		cxi   = np.cos(xi/raddeg)    # cos(xi)
		sxi   = np.sin(xi/raddeg)    # sin(xi)
		cxi2  = cxi * cxi   # cos^2(xi)
		sxi2  = sxi * sxi   # sin^2(xi)
		cxi3  = cxi * cxi2  # cos^3(xi)
		sxi3  = sxi * sxi2  # sin^3(xi)
		cxi4  = cxi2 * cxi2 # cos^4(xi)
		sxi4  = sxi2 * sxi2 # sin^4(xi)
		cxi5  = cxi * cxi4  # cos^5(xi)
		sxi5  = sxi * sxi4  # sin^5(xi)
		cxi6  = cxi2 * cxi4 # cos^6(xi)
		sxi6  = sxi2 * sxi4 # sin^6(xi)
		cxi8  = cxi4 * cxi4 # cos^8(xi)
		sxi8  = sxi4 * sxi4 # sin^8(xi)
		cxi10 = cxi6 * cxi4 # cos^10(xi)
		sxi10 = sxi6 * sxi4 # sin^10(xi)
		cxi12 = cxi6 * cxi6 # cos^12(xi)
		sxi12 = sxi6 * sxi6 # sin^12(xi)  

		psk2 = (1. + 2. * cxi2)/4./sxi2      # for the expansion of S(K)
		psk4 = (1. - 20. * cxi2 - 8. * cxi4)/64./sxi4
		psk6 = -(1. - 18. * cxi2 - 88. * cxi4 - 16. * cxi6)/256./sxi6
		psk8 = -(7. + 8. * cxi2 + 4208. * cxi4 + 5952. * cxi6 + 640. * cxi8)/16384./sxi8
		psk10 = (11. - 130. * cxi2 + 6000. * cxi4 + 35168. * cxi6 + 24704. * cxi8 + 1792. * cxi10)/65536./sxi10
		psk12 = (21. - 204. * cxi2 - 16984. * cxi4 - 414528. * cxi6 - 946560. * cxi8 \
				- 406016. * cxi10  - 21504.* cxi12)/1048576./sxi12	
		

		'''
		Polynomials for the different orders of the solution of the precession angle.
		  	pijpk or pijmk  : i: order of the term in the expansion 
			   				  j: j*zeta in the angular argmuent 
					    pk (mk): +k*alpha (-k*alpha)[p: plus, m: minus] in the angular argument.
		example : p33m4  would mean coefficient in the terms of third order v3 of the trigonometric line 
					with (3*zeta-4*alph	
		'''
		p11p1 = -cxi
		p10p1 = cxi 

		p21p0 = cxi2/2.
		p21p2 = cxi2/2.
		p22p2 = -(1 + 2 * cxi2)/8.
		p20p2 = (1 - 2 * cxi2)/8. 

		p31m1 = cxi/16.
		p31p1 = cxi/16.
		p31p3 = -(cxi - 4. * cxi3)/16.
		p32p1 = -(cxi + 2. * cxi3)/8.
		p32p3 = -(cxi + 2. * cxi3)/8.
		p33p3 = (5. * cxi + 4. * cxi3)/48.
		p30p1 = cxi3/4.
		p30p3 = (cxi - cxi3)/12.  

		p41p0 = -(cxi2 - 4. * cxi4)/16.
		p41m2 = cxi2/96.
		p41p2 = cxi4/8.
		p41p4 = (7. * cxi2 - 12. * cxi4)/96.
		p42p0 = (1. + 4. * cxi2 + 4. * cxi4)/64.
		p42p2 = -(3. - 8. * cxi2 - 6. * cxi4)/48.
		p42p4 = -(1. - 4. * cxi2 - 12. * cxi4)/64.
		p43p2 = -(5. * cxi2 + 4. * cxi4)/32.
		p43p4 = -(5. * cxi2 + 4. * cxi4)/32.
		p44p4 = (3. + 52. * cxi2 + 24 * cxi4)/768.
		p40p2 = (1. - 2. * cxi4)/16.
		p40p4 = (3. - 12. * cxi2 + 8. * cxi4)/256. 	

		#Polynomials for spin phase

		#Secular terms
		qvs2 = 1./2.
		qvs4 = (1. + 3. * cxi2)/8.
		qvs6 = (1. + 2. * cxi2 + 5. * cxi4)/16.
		qvs8 = (5. + 9. * cxi2 + 15. * cxi4 + 35. * cxi6)/128. 
		#  cxi8 corrected to cxi6 14/05/2012 after mail of L. Lindegren
		
		# Fourier terms qixj  i : order (1..5), x = c or s for Cosine or Sine terms, 
		#                     j = Fourier order in the cos or sin
		q1c1 = sxi2
		q2s2 = -sxi2 * cxi/4.
		q3c1 = sxi2 * (1. + 6. * cxi2)/8.
		q3c3 = sxi2 * (1. - 2. * cxi2)/24.
		q4s2 = -sxi2 * cxi3/4.
		q4s4 = -sxi2 * (cxi - cxi3)/32.
		q5c1 = sxi2 * (3. + 12. * cxi2 + 40. * cxi4)/64.
		q5c3 = sxi2 * (9. + 12. * cxi2 - 40. * cxi4)/384.
		q5c5 = sxi2 * (3. - 12. * cxi2 + 8. * cxi4)/640.
		# New to S^7 20/12/2015 to better fit the numerical solution 
		q6s2 = -sxi2 *(-cxi + 2.*cxi3 + 15.*cxi5)/64.
		q6s4 = -sxi2 *( cxi + 2.*cxi3 - 3.*cxi5)/64.
		q6s6 = -sxi2 *( cxi - 2.*cxi3 + cxi5)/192.
		q7c1 = sxi2 * (25. + 90. * cxi2 + 200. * cxi4 + 560. * cxi6)/1024.
		q7c3 = sxi2 * (45. + 90. * cxi2 +  40. * cxi4 - 336. * cxi6)/3072.
		q7c5 = sxi2 * (25. - 30. * cxi2 - 120. * cxi4 + 112. * cxi6)/5120.
		q7c7 = sxi2 * (5.  - 30. * cxi2 +  40. * cxi4 -  16. * cxi6)/7168.
		saa = xi	
		
		parxi = xi, cxi, sxi, cxi2, sxi2, cxi3, sxi3, cxi4, sxi4, cxi5, sxi5, cxi6, sxi6, cxi8,\
				sxi8, cxi10, sxi10, cxi12, sxi12, psk2,\
				psk4, psk6, psk8, psk10, psk12, p11p1, p10p1, p21p0, p21p2, p22p2, p20p2, p31m1,\
				p31p1, p31p3, p32p1, p32p3, p33p3, p30p1, p30p3, p41p0, p41m2, p41p2, p41p4, p42p0,\
				p42p2, p42p4, p43p2, p43p4, p44p4, p40p2, p40p4,\
				qvs2, qvs4, qvs6, qvs8, q1c1, q2s2, q3c1, q3c3, q4s2, q4s4, q5c1, q5c3, q5c5,\
				q6s2, q6s4, q6s6, q7c1, q7c3, q7c5, q7c7, saa
		cputt[0] += time.time()# ctt		
	cputt[1]-=time.time()# ctt
	#Parameters below moved out the precomputed quantities to allow a change of precession rate during the mission. 
	#FM 21 January 2017  & 14 March 2018 for xk

	xk2 = xk * xk
	xk4 = xk2 * xk2
	xk6 = xk2 * xk4
	xk8 = xk4 * xk4
	xk10 = xk8 * xk2
	xk12 = xk8 * xk4
	
	# uniform speed of the z-axis in mean solar motion units
	xs  = xk * sxi * (1. + psk2/xk2 + psk4/xk4 + psk6/xk6 + psk8/xk8 + 1.*psk10/xk10 + 1.*psk12/xk12) 
	xs2 = xs * xs
	xs3 = xs2 * xs
	xs4 = xs2 * xs2
	xs5 = xs3 * xs2
	xs6 = xs3 * xs3
	xs7 = xs4 * xs3
	xs8 = xs4 * xs4
	zaxis_speed = xs # Global variables for the scanning law accessible in the module
	cputt[1]+=time.time()# ctt
	#  Second group : depending on the intial condition of the scans. This can be changed during a run.
	cputt[9]-=time.time()# ctt
	cketaz, sketaz =  trigmul(etaz/raddeg, imax)
	cputt[9]+=time.time()# ctt
	''' ------------------------------------------------------------------------------------------------------
	motion of the spin axis on the sky - Euler anglespsi and teta
	sun longitude in degrees. Conventional sun for the scanning law.
	'''
	cputt[2]-=time.time()# ctt
	xlsun_ref = sun_long(xjd_ref, sun_ra)
	cputt[2]+=time.time()# ctt

	cputt[3]-=time.time()# ctt
	xjd = date1 + date2
	xlsun = sun_long(xjd, sun_ra)
	cputt[3]+=time.time()# ctt
	#------------------------------------------------------------------------------------------------------

	# ctt
	'''------------------------------------------------------------------------------------------------------
	eta: precession angle on the sun-centered cone at xmjd in deg
		 Denoted v in FM-037, FM-059
	'''
	arg = (xlsun - xlsun_ref) * xk
	cputt[9]-=time.time()
	ckarg, skarg  = trigmul(arg/raddeg, 4)
	cputt[9]+=time.time()
	cputt[4]-=time.time()
	# zero order
	eta0 = arg + etaz # this is zeta + alpha in the TN

	# first order
	cc, ss = trigsum(ckarg, skarg, cketaz, sketaz)
	eta1 = (p11p1 * cc[1,1] + p10p1 * cc[0,1])/xs

	#eta1 = -cos(xi)/xs (cc - cetaz[1])
	# second order
	eta2 = p21p0 * ss[1,0]
	eta2 = eta2 + p21p2 * ss[1,2]
	eta2 = eta2 + p22p2 * ss[2,2]
	eta2 = eta2 + p20p2 * sketaz[2]
	eta2 = eta2/xs2

	# third order

	eta3 = p31m1 * cc[1,-1]
	eta3 = eta3 + p31p1 * cc[1,1]
	eta3 = eta3 + p31p3 * cc[1,3]
	eta3 = eta3 + p32p1 * cc[2,1]
	eta3 = eta3 + p32p3 * cc[2,3]
	eta3 = eta3 + p33p3 * cc[3,3]
	eta3 = eta3 + p30p1 * cc[0,1]
	eta3 = eta3 + p30p3 * cc[0,3]

	eta3 = eta3/xs3
	# fourth order
	eta4 = p41p0 * skarg[1]
	eta4 = eta4 + p41m2 * ss[1,-2]
	eta4 = eta4 + p41p2 * ss[1,2]
	eta4 = eta4 + p41p4 * ss[1,4]
	eta4 = eta4 + p42p0 * ss[2,0]
	eta4 = eta4 + p42p2 * ss[2,2]
	eta4 = eta4 + p42p4 * ss[2,4]
	eta4 = eta4 + p43p2 * ss[3,2]
	eta4 = eta4 + p43p4 * ss[3,4]
	eta4 = eta4 + p44p4 * ss[4,4]
	eta4 = eta4 + p40p2 * ss[0,2]
	eta4 = eta4 + p40p4 * ss[0,4]
	eta4 = eta4/xs4
	# full solution in deg


	eta = eta0 + raddeg * (eta1 + eta2 + eta3 + eta4) 
	cputt[4]+=time.time()# ctt
	
	cputt[9]-=time.time()
	cketa, sketa =  trigmul(eta/raddeg, imax)
	cputt[9]+=time.time()


	#--------------------------------------------------------------	
	cputt[5]-=time.time()# ctt

	''' --------------------------------------------------------------
	Patch for numerical complements to numerical scanning law - 
	Numerical values for SL_SMOOTH, GAREQ, Reverse_SL, post_reverse only - 17/12/2015, May 2018, May 2019  
	Numerical fitting from differences between Numerical and analytical SL
	Amplitudes in arcsec
	'''
	domegaz = 0e0
	deta = 0e0
	
	if(abs(etaz-52.6600) <1e-2):   # 52.66 : initial Spin Precession Gareq

		camp[1] =  4.46195e0
		camp[2] = -9.05210e0
		camp[3] = -3.61607e0
		camp[4] =  9.82721e0
		camp[5] = -1.63379e0
		camp[6] = -0.63960e0
		camp[7] =  0.12777e0
		camp[8] =  0.0041739e0 

		samp[1] =  4.88724e0
		samp[2] = -0.736780e0
		samp[3] = -5.24563e0
		samp[4] =  5.79476e0
		samp[5] = -3.65478e0
		samp[6] =  0.781132e0
		samp[7] =  0.0638473e0
		samp[8] = -0.010442e0 

		deta = -1.19784e0 - 7.17091e-6*((date1-xjd_ref) +date2)
		for ic in range(1, 9):
			deta = deta + camp[ic]*cketa[ic] + samp[ic]*sketa[ic]

		domegaz = 4.7702e-9  # Drift between numerical and analytical in deg/hour
	elif (abs(etaz-270.65) <1e-2): # 270.65: initial Spin Precession for the reverse SL (May 27, 2019)
		camp = np.zeros(9) # in arcsec
		camp[1] =  -2.44981e0
		camp[2] =   0.19465e0
		camp[3] =  -7.69931e0
		camp[4] =  -0.20277e0
		camp[5] =   2.81329e0
		camp[6] =   0.01239e0
		camp[7] =  -0.14181e0
		camp[8] =  -3.14974e-5 

		samp = np.zeros(9) #in arcsec
		samp[1] =  0.11584e0
		samp[2] = -3.24763e0
		samp[3] = -0.04175e0
		samp[4] =  0.48111e0
		samp[5] = -0.07184e0
		samp[6] =  0.93939e0
		samp[7] =  0.00121e0
		samp[8] = -0.00989e0 

		deta = 0.061406e0 +  7.17914e-6*((date1-xjd_ref) +date2)
		for ic in range(1, 9):
			deta = deta + camp[ic]*cketa[ic] + samp[ic]*sketa[ic]

		domegaz = -4.7702e-9   # Drift between numerical and analytical in deg/hour

	elif (abs(etaz-204.64) <1e-2): # 204.64: initial Spin Precession for the extension after reverse SL (May 27, 2019)
		camp = np.zeros(9) # in arcsec
		camp[1] =  10.20026e0
		camp[2] =  10.81124e0
		camp[3] = -16.65216e0
		camp[4] = -12.85838e0 
		camp[5] =  -0.21342e0
		camp[6] =   0.91408e0
		camp[7] =   0.10705e0
		camp[8] =  -0.00877e0

		samp = np.zeros(9) # in arcsec
		samp[1] =  -0.44130e0
		samp[2] =   3.59759e0
		samp[3] =  13.39774e0
		samp[4] =  12.12510e0  
		samp[5] =   5.13206e0
		samp[6] =   0.58803e0
		samp[7] =  -0.09332e0
		samp[8] =  -0.00861e0

		deta = -1.59021e0 - 7.14559e-6*((date1-xjd_ref) +date2)
		for ic in range(1, 9):
			deta = deta + camp[ic]*cketa[ic] + samp[ic]*sketa[ic]
			
		domegaz = 4.7702e-9  # Drift between numerical and analytical in deg/hour

	elif(abs(etaz- 290.745)<1e-2): # 290.745D0 Initial precession phase for smooth NSL        
		camp = np.zeros(9) # in arcsec
		camp[1] =  2.8769
		camp[2] = -5.3119
		camp[3] =  3.7565
		camp[4] =  6.1317 
		camp[5] = -2.4095
		camp[6] = -0.3813
		camp[7] =  0.1374
		camp[8] =   0.0017

		samp = np.zeros(9) # in arcsec
		samp[1] = -2.4317
		samp[2] = -2.4873
		samp[3] = -1.9070
		samp[4] =  2.3111
		samp[5] = -2.2039
		samp[6] =  0.8856
		samp[7] =  0.0372
		samp[8] = -0.0103

		deta = -3.8301174  -7.10705e-6*((date1-xjd_ref) + date2) # fit over the first 500 days only
		for ic in range(1, 9):
			deta = deta + camp[ic]*cketa[ic] + samp[ic]*sketa[ic]
		      

		domegaz = 4.7682e-9   # Drift between numerical and analytical in deg/hour         

	eta = eta + deta*asdeg # add the numerical complements - comment the statement to disable the numerical complements
	cputt[5]+=time.time()# ctt

	#--------------------------------------------------------------

	#--------------------------------------------------------------
	cputt[6]-=time.time()# ctt
	ceta = np.cos(eta/raddeg)
	seta = np.sin(eta/raddeg)
	#
	#Auxiliary quantities
	#
	x = cxi
	y = sxi * ceta
	z = sxi * seta

	psi = np.arctan2(x, -y) # first  Euler angle
	teta = np.arccos(z)     # second Euler angle
	cputt[6]+=time.time()# ctt
	cputt[7]-=time.time()# ctt

	stet = np.sin(teta)
	ctet = np.cos(teta)
	cpsi = np.cos(psi)
	spsi = np.sin(psi)
	cputt[7]+=time.time()# ctt

	cputt[6]-=time.time()# ctt
	psi  = psi * raddeg  # node/sun
	teta = teta * raddeg # inclination
	cputt[6]+=time.time()# ctt

	#    Time derivatives in deg/hour (same as "/s)
	cputt[7]-=time.time()# ctt
	xlamdap = dsl_dt(xjd, sun_ra)/dayh # angular motion of the sun in longitude in deg/h
	# Equation for etap changed on 30 Oct 2019 - Before xs was assumed to be >0 
	# and the result was wrong for reverse SL (wrong also for psip, phip, tetap)

	# angular speed in "/s or deg/h for the precession angle from the defining equation
	etap = 1./sxi * xlamdap * (xs*np.sqrt((1. - ceta*ceta/xs2)) + cxi * seta) 
	psip = -spsi * ctet/stet * etap + xlamdap # includes the solar motion (INERTIAL psip and not wrt solar node)
	tetap = cpsi * etap
	#
	#     projected sun  to node in the scanning circle :: heliotropic longitude of the node
	#     two formula with identical results	

	# omenode stands for 'omega of node'"
	omenode = np.arctan2(spsi * ctet, cpsi) * raddeg
	#
	#     rotation about the spin axis   -  heliotropic value omega and euler angle phi 
	#
	s2eta = 2. * ceta * seta
	c2eta = 2. * ceta * ceta - 1.
	s3eta = s2eta * ceta + seta * c2eta
	c3eta = c2eta * ceta - s2eta * seta
	s4eta = 2. * c2eta * s2eta
	c4eta = 2. * c2eta * c2eta - 1.
	s5eta = s4eta * ceta + seta * c4eta
	c5eta = c4eta * ceta - s4eta * seta
	s6eta = s4eta * c2eta + s2eta * c4eta
	c6eta = c4eta * c2eta - s4eta * s2eta
	s7eta = s6eta * ceta + seta * c6eta
	c7eta = c6eta * ceta - s6eta * seta
	s8eta = 2. * s4eta * c4eta
	c8eta = c4eta * c4eta - s4eta * s4eta
	#
	#Recoding with order 1/S^5 03/04/2012
	#
	omega = phiz  \
		+(date1 - xjd_ref) * ((omegaz +domegaz) * dayh) + date2 * ((omegaz +domegaz) * dayh) \
		+sxi2 * cxi * (qvs2/xs2 + qvs4/xs4 + qvs6/xs6 + qvs8/xs8)*(eta - etaz) \
		-cxi * (eta - etaz) \
		+raddeg * (\
		+1./xs * q1c1 * (ceta - cketaz[1]) \
		+1./xs2 * q2s2 * (s2eta - sketaz[2])\
		+1./xs3 * (q3c1 * (ceta - cketaz[1]) + q3c3 * (c3eta - cketaz[3])) \
		+1./xs4 * (q4s2 * (s2eta - sketaz[2]) + q4s4 * (s4eta - sketaz[4]))  \
		+1./xs5 * (q5c1 * (ceta - cketaz[1]) + q5c3 * (c3eta - cketaz[3]) + q5c5 * (c5eta - cketaz[5])) \
		+1./xs6 * (q6s2 * (s2eta - sketaz[2]) + q6s4 * (s4eta - sketaz[4])+ q6s6 * (s6eta - sketaz[6]))  \
		+1./xs7 * (q7c1 * (ceta - cketaz[1]) + q7c3 * (c3eta - cketaz[3]) + q7c5 * (c5eta - cketaz[5])\
		+ q7c7 * (c7eta - cketaz[7])) \
	) # end of raddeg

	if  xk  > 0:  
		# direct SL  (amplitudes in arcsec)
		omega = omega +   (  0.0298*(ceta-cketaz[1])   +  (-0.0172)*(s2eta- sketaz[2]) \
			+ 0.0052*(c3eta- cketaz[3]) + (-0.01415)*(s4eta-sketaz[4]) + (-0.00204)*(s6eta- sketaz[6]))*asdeg
	else:               
		#reverse SL
		omega = omega +   ((-0.0298)*(ceta-cketaz[1])  +  (-0.0172)*(s2eta- sketaz[2]) \
			+ (-0.0052)*(c3eta- cketaz[3]) + (-0.01415)*(s4eta-sketaz[4]) + (-0.00204)*(s6eta- sketaz[6]))*asdeg
	omega = omega % ro  # commented for thermistor investigation  Nov 2015 - Should be unchecked
	phi = (omega - omenode) % ro # phi is measured from the node

	#domega/dt in deg/h - exact value from the defining equation
	omegap = omegaz +domegaz - xlamdap * sxi * seta - etap * cxi 
	phip = omegaz - xlamdap * sxi * seta + etap * (cxi * sxi2 * seta * seta)/(1. - sxi2 * seta * seta)

	cputt[7]+=time.time()# ctt
	
	cputt[8]-=time.time()# ctt
	'''------------------------------------------------------------------------------------------------------
	calculate the Coordinates of the two fied of views
	------------------------------------------------------------------------------------------------------'''
	sxlsun,cxlsun = np.sin(np.pi*xlsun/180), np.cos(np.pi*xlsun/180)
	somega,comega = np.sin(omega/raddeg), np.cos(omega/raddeg)
	cba,sba = csba
	#Cartesian coordinates z axis:
	z_x = sxlsun * cxi + ceta * cxlsun * sxi
	z_y = cxlsun * cxi - ceta * sxlsun * sxi 
	z_z	= seta * sxi
	z = np.array([z_x,z_y,z_z])
	


	#Cartesian coordinate preceding FOV
	pfov_x =  (comega * cba - somega * sba) * (sxlsun * sxi - ceta * cxlsun * cxi) \
				- (somega * cba + comega * sba) * seta * cxlsun
	pfov_y = (comega * cba - somega * sba) * (cxlsun * sxi + ceta * sxlsun * cxi) \
				+ (somega * cba + comega * sba) * seta * sxlsun
	pfov_z = -(comega * cba - somega * sba) * (seta * cxi) + (somega * cba + comega * sba) * ceta	
	#Cartesian coordinate following FOV
	ffov_x =  (comega * cba + somega * sba) * (sxlsun * sxi - ceta * cxlsun * cxi) \
				- (somega * cba - comega * sba) * seta * cxlsun
	ffov_y = (comega * cba + somega * sba) * (cxlsun * sxi + ceta * sxlsun * cxi) \
				+ (somega * cba - comega * sba) * seta * sxlsun
	ffov_z = -(comega * cba + somega * sba) * (seta * cxi) + (somega * cba - comega * sba) * ceta

	#calculate position angle of scan-direction
	pfov_dir_x = -(somega * cba + comega * sba) * (sxlsun * sxi - ceta * cxlsun * cxi) \
				- (comega * cba - somega * sba) * seta * cxlsun
	pfov_dir_y = -(somega * cba + comega * sba) * (cxlsun * sxi + ceta * sxlsun * cxi) \
				+ (comega * cba - somega * sba) * seta * sxlsun
	
	pfov_dir_z = (somega * cba + comega * sba) * (seta * cxi) + (comega * cba - somega * sba) * ceta
	
	ffov_dir_x =  (- somega * cba + comega * sba) * (sxlsun * sxi - ceta * cxlsun * cxi) \
				- (comega * cba +somega * sba) * seta * cxlsun
	ffov_dir_y = (- somega * cba + comega * sba) * (cxlsun * sxi + ceta * sxlsun * cxi) \
				+ (comega * cba +somega * sba) * seta * sxlsun
	
	ffov_dir_z = (somega * cba - comega * sba) * (seta * cxi) + (comega * cba + somega * sba) * ceta

	#Translate Coordinates

	ffov = np.array([ffov_x,ffov_y,ffov_z])
	pfov = np.array([pfov_x,pfov_y,pfov_z])
	ffov_dir = np.array([ffov_dir_x,ffov_dir_y,ffov_dir_z])
	pfov_dir = np.array([pfov_dir_x,pfov_dir_y,pfov_dir_z])
		#------------------------------------------------------------------------------------------------------
	cputt[8]+=time.time() # ctt
	
	#return xlsun,eta, psi, teta, omega, phi, etap, psip, tetap, omegap ,phip, omenode, z,pfov,ffov, parxi, tt
	return {'xlsun':xlsun, 'eta':eta, 'psi':psi, 'teta':teta, 'omega':omega, 'phi':phi, 'etap':etap, 'psip':psip,\
			'tetap':tetap, 'omegap':omegap ,'phip':phip, 'omenode':omenode, 'date':date1+date2, \
			'z':z, 'pfov':pfov, 'pfov_dir':pfov_dir, 'ffov':ffov, 'ffov_dir':ffov_dir,\
			'parxi':parxi, 'cputt':cputt}

def n_row(z,pos, hight = 0.60):
	'''------------------------------------------------------------------------------------------------------
	Determine the ccd Row
	------------------------------------------'''
	zdotpos = z.dot(pos)
	raddeg  = 180 / np.pi 	  						# radians to degree

	if  zdotpos>np.sin((hight/2-hight/7)/raddeg): return 1
	elif  zdotpos>np.sin((hight/2-2*hight/7)/raddeg): return 2
	elif  zdotpos>np.sin((hight/2-3*hight/7)/raddeg): return 3
	elif  zdotpos>np.sin((hight/2-4*hight/7)/raddeg): return 4
	elif  zdotpos>np.sin((hight/2-5*hight/7)/raddeg): return 5
	elif  zdotpos>np.sin((hight/2-6*hight/7)/raddeg): return 6
	else: return 7

def scandir(z,fov, fov_dir):
	'''------------------------------------------------------------------------------------------------------
	Calculate the position angle in radian of the scan direction from cartesian vectors of the gaia aglinment
	------------------------------------------'''
	if z[2]>0:
		scan_dir = np.arccos(fov_dir[2]/np.sqrt(1-fov[2]*fov[2]))
	elif z[2]==0:
		if fov_dir[2] >0:  scan_dir = 0
		else: scan_dir = np.pi
	else: #z[2]<0
		scan_dir =  2*np.pi-np.arccos(fov_dir[2]/np.sqrt(1-fov[2]*fov[2]))
	return scan_dir

def vary_eta(starlist, tmin, tmax, etaz=0, eps_deg=0.5,eps_ac = 0.30, grid_s=120, grid_d=2, parxi=None, sun_ra = None, **kwargs):
	print(etaz)
	'''------------------------------------------------------------------------------------------------------
	Calculate the algliment of Gaia for a given initial eta.
	Uses a starlist to finde important epochs within a sparse grid,
	Around thes epochs a dens grid is used.	
	---------------------------------------------------------------------------------------------------------

		Input
			starlist  : coordinates of the test stars to reduce the computationtime
			tmin; tmax:	time range for calculation				 
			etaz 	  :	initial phase of the precession motion  in deg at xjd_ref (see scaninglaw)
			eps_deg   : accuracy limit for the z axis alignment(°)
			eps_ac 	  : limit for the x axis alignment(°) (half hight of the fov window )
			grid_s	  :	grid for the sparsely calculations (min)
			grid_d    : grid for the dens calculation (min)
			parxi	  :	precalculated parameters for xi 						  (see scaninglaw)	
			sun_ra 	  : precalculated longitudes of the sun 				      (see scaninglaw)

	---------------------------------------------------------------------------------------------------------
		Output
			alignment_dense : list of gaia alginment in a dense grid at important epochs 
			cputt           : compunting time tracker
			cputt_sclaw     : compunting time tracker of the scaninglaw subroutine
	------------------------------------------------------------------------------------------------------'''
	#compunting-time-tracker
	cputt = np.zeros(6)
	#------------------------------------------------------------------------------------------------------
	#constants
	omegaz = 59.960499070315						# nominal scanning speed in deg/hour (= "/s  ~ 60.0)
	minday  = 1 / 24 / 60 							# minuts in unit days 
	raddeg  = 180 / np.pi 	  						# radians to degree
	eps2    = np.sin(eps_deg / raddeg)**2 			# limit in cartesian coordinates squared for first estimate
	eps_al	= np.cos(omegaz/30 / raddeg)	# limit in cartesian coordinates squared in along-scan dir.
	eps_ac2	= np.sin(eps_ac / raddeg)**2			# limit in cartesian coordinates squared in cross-scan dir.
	#------------------------------------------------------------------------------------------------------

	#------------------------------------------------------------------------------------------------------
	#translate spherical coordinates to cartesian coordinates
	#		cos(delta) sin(alpha)
	#	x = cos(delta) cos(alpha)
	#		sin(delta) 
	#
	#cputt[0]-=time.time()
	starlist_name = starlist['source_id']
	starlist_cart = np.array([np.sin(starlist['ra']/raddeg)*np.cos(starlist['dec']/raddeg), \
							  np.cos(starlist['ra']/raddeg)*np.cos(starlist['dec']/raddeg), \
							  np.sin(starlist['dec']/raddeg)]).T
	#cputt[0]+=time.time()
	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	# create sparse time grid
	#cputt[1]-=time.time()
	xjd_ref=Time(tmin,format='jyear').jd # reference epoch
	t1 = np.arange(Time(tmin,format='jyear').jd,Time(tmax,format='jyear').jd,minday*grid_s) # time_grid
	#cputt[1]+=time.time()
	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	#calculate alignment of gaia 
	cputt[2] -= time.time()
	alignment = [scaninglaw(date1=xx, xjd_ref=xjd_ref, etaz=etaz, parxi = parxi, sun_ra = sun_ra) for xx in t1]
	cputt[2] += time.time()
	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	#finde timeranges when one source is perpendiculare to the z-axis
	cputt[2] -= time.time()
	alignment_filterd = list(filter(lambda z: (starlist_cart.dot(z['z'])**2 < eps2).any(), alignment))
	cputt[2] += time.time()
	#------------------------------------------------------------------------------------------------------

	#------------------------------------------------------------------------------------------------------
	# create dense time grid
	cputt[3] -= time.time()
	t2 = np.array([x['date'] for x in alignment_filterd]) 
	dt = np.arange(-grid_s/2.+grid_d/2.,grid_s/2.+grid_d/2.,grid_d)
	t2 = np.repeat(t2.reshape(-1,1),len(dt), axis = 1)
	dt =  np.repeat(dt.reshape(1,-1),len(alignment_filterd), axis= 0)
	t2 = (t2+minday*dt).reshape(-1)
	cputt[3] += time.time()
	#------------------------------------------------------------------------------------------------------

	#------------------------------------------------------------------------------------------------------
	#calculate alignment of gaia 
	cputt[4]-= time.time()
	alignment_dense = [scaninglaw(date1=xx, xjd_ref=xjd_ref, etaz=etaz, parxi = parxi, sun_ra = sun_ra) for xx in t2]
	cputt[4]+= time.time()
	#------------------------------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------------------------
	#find observation epochs

	cputt[5] -= time.time()
	tab = []
	for star_pos,star_name in zip(starlist_cart,starlist_name):
		#filter epochs
		pfov_filterd = np.array(list(filter(lambda z: (star_pos.dot(z['pfov']) > eps_al) and (star_pos.dot(z['z'])**2 < eps_ac2), alignment_dense)))
		ffov_filterd = np.array(list(filter(lambda z: (star_pos.dot(z['ffov']) > eps_al) and (star_pos.dot(z['z'])**2 < eps_ac2), alignment_dense)))		
		#correct observation time
		t_pfov       = np.array([x['date']+np.arcsin(star_pos.dot(x['pfov_dir']))*raddeg/omegaz/24 for x in pfov_filterd])
		t_ffov       = np.array([x['date']+np.arcsin(star_pos.dot(x['ffov_dir']))*raddeg/omegaz/24 for x in ffov_filterd])
		#remove doubled observations
		pfov_filterd = pfov_filterd[np.where(np.append(t_pfov[1:],9999999999)>t_pfov+0.125)]
		t_pfov       = t_pfov[np.where(np.append(t_pfov[1:],9999999999)>t_pfov+0.125)]
		ffov_filterd = ffov_filterd[np.where(np.append(t_ffov[1:],9999999999)>t_ffov+0.125)]
		t_ffov       = t_ffov[np.where(np.append(t_ffov[1:],9999999999)>t_ffov+0.125)]
		#calculate scandir and ccd_row
		if len(pfov_filterd)>0:
			sdn_pfov,nrow_pfov = np.array([(scandir(x['z'],x['pfov'],x['pfov_dir']),n_row(x['z'],star_pos)) for x in pfov_filterd]).T
		else:
			sdn_pfov = []
			nrow_pfov = [] 

		if len(ffov_filterd)>0:
			sdn_ffov,nrow_ffov = np.array([(scandir(x['z'],x['ffov'],x['ffov_dir']),n_row(x['z'],star_pos)) for x in ffov_filterd]).T
		else:
			sdn_ffov = []
			nrow_ffov = []
		#------------------------------------------------------------------------------------------------------
		#Create table 	
		tab_i = Table()
		tab_i['Target'] = [star_name]*(len(sdn_pfov)+len(sdn_ffov))
		tab_i['scanAngle[rad]'] = np.append(sdn_pfov,sdn_ffov)
		tab_i['CcdRow[1-7]'] = np.append(nrow_pfov,nrow_ffov)
		tab_i['ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]'] = np.append(t_pfov,t_ffov)
		tab_i['eta[degree]'] = np.ones(len(tab_i))*etaz
		tab.append(tab_i)
	tab = vstack(tab)
	cputt[5] += time.time()

	cputt_sclaw1 = np.array([a['cputt'] for a in alignment])
	cputt_sclaw2 = np.array([a['cputt'] for a in alignment_dense])
	cputt_sclaw =  np.sum(cputt_sclaw1, axis=0) + np.sum(cputt_sclaw2, axis=0)
	return tab, cputt, cputt_sclaw

def create_event_list(dist_lim=350, tmin=2020.5, tmax=2023.5, GL_low=7, dG_up=5, shift_lim = 0.25, shift_lim_max=[0.5,0.1]):
	'''------------------------------------------------------------------------------------------------------
	Filter the eventlist for dozen events for optimisation
	------------------------------------------'''
	
	tab = Table.read(path + 'InputTable/events_within_gaia.fits', format = 'fits')
	tab = tab[np.where(tab['L2_date']>tmin)]
	tab = tab[np.where(tab['L2_date']<tmax)]
	tab = tab[np.where(tab['phot_g_mean_mag']> GL_low)]
	tab = tab[np.where(tab['ob_phot_g_mean_mag']-tab['phot_g_mean_mag'] <dG_up)]
	z = 10 ** (0.4 * (np.maximum(tab['ob_phot_g_mean_mag'], 14) - 15))
	sig_fov =((-1.631 + 680.766 * z + 32.732 * z**2)**0.5/1000 *7.75+0.1)
	shift = np.maximum(tab['L2_dist'],dist_lim)*tab['thetaE']*tab['thetaE']\
			/ (np.maximum(tab['L2_dist'],dist_lim)*np.maximum(tab['L2_dist'],dist_lim)+2*tab['thetaE']*tab['thetaE'])
	tab1 = tab[np.where((shift > shift_lim) | ((tab['L2_shift_plus'] > shift_lim_max[0]) & (shift > shift_lim_max[1])))]
	tab2 = tab.copy()
	tab2['shift_over_sig'] = shift/sig_fov*3
	tab2.sort('shift_over_sig')
	tab2.reverse()
	for sid in np.unique(tab2['source_id']):
		where = np.where(tab2['source_id'] == sid)[0]
		if len(where)>1:
			tab2.remove_rows(where[1:])
			print(len(np.where(tab2['source_id'] == sid)[0]))

	print(len(tab1))
	tab2 = tab2[np.where(tab2['shift_over_sig']>1)]
	for i in range(len(tab2)-1,-1,-1):
		if tab2['source_id'][i] in tab1['source_id']:
			tab2.remove_rows([i])

	print(len(tab2))
	tab1.write(path + 'InputTable/events_eta.fits', format = 'fits', overwrite = True)
	tab2.write(path + 'InputTable/events_eta_2.fits', format = 'fits', overwrite = True)

def vary_eta_paralle(i,parts,starlist,keywords={}):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	d_eta = keywords.get('d_eta', 90)
	eta = np.arange(i*360/parts,(i+1)*360/parts,d_eta)
	tab,cputt_vary,cputt_sclaw = \
			zip(*[vary_eta(starlist, 2020.5, 2024,etaz = e, **keywords) for e in eta])
	tab = vstack(tab)
	cputt_vary = np.array(cputt_vary)
	cputt_vary = np.sum(cputt_vary, axis = 0)
	cputt_sclaw = np.array(cputt_sclaw)
	cputt_sclaw = np.sum(cputt_sclaw, axis = 0)
	return tab, cputt_vary,cputt_sclaw
def main(test = False, tab = 1, part=[0,0]):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	if 'jonas' in path: 
		test = True 
		print(test)
	'''------------------------------------------------------------------------------------------------------
	Calculate the t
	------------------------------------------------------------------------------------------------------'''
	#compunting-time-tracker
	cputt = np.zeros(5)
	#------------------------------------------------------------------------------------------------------'''
	#constants
	d_eta   = 0.05							# steps to vary eta (degree)
	if test: d_eta = 90 
	minday  = 1 / 24 / 60 					# minuts in unit days 
	hday    = 1 / 24.						# hours in unit days
	if test: hday = hday*200 
	raddeg  = 180 / np.pi 	  				# radians to degree
	grid_s  = 120							# grid densitiy sparsely (min)
	grid_d  = 2 							# grid densitiy dens (min)
	if test: grid_d = grid_d*20 
	eps_deg = 0.5							# limit for the z axis in degree()
	eps2    = np.sin(eps_deg / raddeg)**2 	# limit in cartesian coordinates squared 
	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	#loadstarlist
	#cputt[0]-=time.time()
	if tab <= 1:
		tabstring = ''
		starlist = Table.read(path+'InputTable/events_eta.fits', format = 'fits')
	elif tab <= 2:
		tabstring = '2'
		starlist = Table.read(path+'InputTable/events_eta_2.fits', format = 'fits')
	else:
		tabstring = '3'
		starlist = Table.read(path+'InputTable/events_eta_3.fits', format = 'fits')

	print(len(starlist))

	#cputt[0]+=time.time()
	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	#calculate sun_orbit and parameters once to pass to subroutines
	cputt[1]-=time.time()
	t_ra    = np.arange(Time(2020, format='jyear').jd,Time(2025, format='jyear').jd,hday)	#1 hour grid
	xjd_ref = Time(2020.5, format='jyear').jd 												#reference epoch
	xlsun   = getSun_long_vec(t_ra)														#longitudes of the sun
	sun_ra  = (t_ra,xlsun)														
	parxi   = scaninglaw(date1=xjd_ref, xjd_ref=xjd_ref, sun_ra = sun_ra)['parxi'] 			#calculate parameters
	cputt[1]+=time.time()
	#------------------------------------------------------------------------------------------------------

	#------------------------------------------------------------------------------------------------------
	# calulate alignment while vary eta
	cputt[2]-=time.time()
	if part[0] == -1:
		keywords= dict(grid_d = 2,d_eta = d_eta, parxi=parxi,sun_ra=sun_ra)
		from  joblib import Parallel, delayed
		tab_1 = Parallel(n_jobs=part[1])(delayed(vary_eta_paralle)(i,part[1],starlist, keywords) for i in range(part[1]))
		cputt_vary = cputt_sclaw = np.zeros((1,1))
		tab = [i[0] for i in tab_1]
		cputt_vary =[i[1] for i in tab_1]
		cputt_sclaw=[i[2] for i in tab_1]
		print(len(tab))
	else:
		if part[1] == 0:
			eta = np.arange(0,360,d_eta) 
		else:
			eta = np.arange(part[0]*360/part[1],(part[0]+1)*360/part[1],d_eta)   #grid on eta
		print(len(eta))
		tab,cputt_vary,cputt_sclaw = \
			zip(*[vary_eta(starlist, 2020.5, 2024, grid_d = 2, etaz=e, parxi=parxi,sun_ra=sun_ra) for e in eta])
	cputt[2]+=time.time()
	tab = vstack(tab)
	print (len(tab))



	#------------------------------------------------------------------------------------------------------

	#------------------------------------------------------------------------------------------------------
	#save results 
	if test:tab.write(path+'InputTable/tab'+tabstring+'_for_vary_eta_part_test.csv', format = 'csv',overwrite=True)
	elif part[0] == -1:	tab.write(path+'InputTable/tab'+tabstring+'_for_vary_eta.csv', format = 'csv',overwrite=True)
	elif part[1] == [0]: 	tab.write(path+'InputTable/tab'+tabstring+'_for_vary_eta.csv', format = 'csv',overwrite=True)
	else:	tab.write(path+'InputTable/tab'+tabstring+'_for_vary_eta_part'+str(part[0])+'.csv', format = 'csv',overwrite=True)


	#contains id, time, theta, ccd row, eta

	#------------------------------------------------------------------------------------------------------
	
	#------------------------------------------------------------------------------------------------------
	# print computingtime trackers:
	print(' main:', cputt)
	cputt_vary = np.array(cputt_vary)
	print(' vary:', np.sum(cputt_vary, axis = 0))
	cputt_sclaw = np.array(cputt_sclaw)
	#print('sclaw:', np.sum(cputt_sclaw, axis = 0))

def analysis(results = None):
	'''------------------------------------------------------------
		Description: 
	---------------------------------------------------------------
		Input:
	---------------------------------------------------------------
		Output:
	------------------------------------------------------------'''
	import MLG
	import matplotlib.pyplot as plt
	from scipy.signal import savgol_filter
	outlist = np.array(Table.read(path + 'out_list', format = 'csv')['source_ID'], str)
	if results is None:
		res_tab = []

		for string in ['MC_single_vary_eta_17-02-2020_1440_10_100_1.pkl',]:#, 'MC_single_vary_eta_17-01-2020_25920_10_500_1.pkl']:
			data = MLG.Simulation.MonteCarlo.loadPkl(string)
			results = MLG.Analysis.statistics(data, string = 'vary_eta')
			res_tab.append(results)
	else: res_tab = results
	delta_eta	 = 0.05
	for i in [2,1]:
		fig = plt.figure(figsize= [15,7])
		ax = plt.subplot(111)
		for tab in range(len(res_tab)):
			ids = np.array(res_tab[tab][1])

			res = np.array([k[0][2:] for k in res_tab[tab][0]]) 
			source_ids = np.unique(ids)
			for si in source_ids:
				which = np.where(ids == si)
				res_si = res[which]
				mass = res_si[:,0]
				plotx = np.array([res_tab[tab][2][j][0].getEta() for j in which[0]])
				if str(si) == '239070631455336064': st = ''
				elif res_tab[tab][2][which[0][0]][1].getPx() == 0:
					st = '!'
					continue
				else : st = ''
				if str(si) == '239070631455336064': pass
				elif si in  outlist: 
					print(si)
					continue
				if res_tab[tab][2][which[0][0]][0].getMass() == 0.65: co = 'red'
				else: co = 'blue' 
				ploty = mass/res_si[:,5]
				#plotx = np.arange(0,360,delta_eta)
				plotx = plotx[np.where(ploty>0)]
				ploty = ploty[np.where(ploty>0)]
				if len(ploty)>51:
					ploty_h= savgol_filter(ploty, 51, 3)
					if i == 1:
						if ((max(ploty)> 4) & (min(ploty)> 1)) : plt.plot(plotx,ploty_h,color = co, label = [si+st,str(mass[0])[:4]])		
					if i == 2:
						if ((max(ploty) > 4) & (min(ploty)> 1)): plt.plot(plotx,ploty,'+',color = co, label =[si+st,str(mass[0])[:4]])
			
		plt.ylabel('M / sig_M')
		plt.xlabel('ny [°]')
		box = ax.get_position()
		plt.plot([0,360],[0.43/0.13,0.43/0.13],'--', color = 'blue',linewidth = 4, label= 'current scanning law')
		plt.plot([0,360],[0.65/0.13,0.65/0.13],'--', color = 'red',linewidth = 4, label= 'current scanning law')
		plt.plot([0,360],[0.65/0.077,0.65/0.077],'k--',linewidth = 4, label= 'result for 4970215770740383616, (T_ca = 2018.098)')
		plt.plot([0,360],[0.65/(0.077+0.007),0.65/(0.077+0.007)],'k--',linewidth = 2)
		plt.plot([0,360],[0.65/(0.077-0.006),0.65/(0.077-0.006)],'k--', linewidth = 2) 
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
		ax.legend(fontsize = 12, loc='center left', bbox_to_anchor=(1, 0.5))

	plt.show()
	return res_tab
def reshape():
	f = open(path +'InputTable/tab3_for_vary_eta_old.csv', 'r')
	g = open(path + 'InputTable/tab3_for_vary_eta.csv', 'w')
	c = 0
	for line in f:
		ll = line
		if c >0:
			ll = ll.split('.')
			if len(ll[-1])> 5: 
				dd = ll[-1]
				if dd[3] == '9':
					if dd[1] == '4':
						ll[-1] = dd[0]+'5\n'
					elif dd[1] == '9': 
						cc = (int(dd[0])+1) % 10
						if cc == 0:
							bb = ll[-2].split(',')
							bb[-1] = str(int(bb[-1])+1)
							ll[-2] = ','.join(bb)
							ll[-1] = '00\n'
						else:
							ll[-1] = str(cc)+'0\n'
					else: print('error2', dd[2])
				elif int(dd[3]) == 0: ll[-1] = dd[:2]+'\n'
				else: print('error3', dd[3])

			ll='.'.join(ll)
		g.write(ll)
		c+=1
	g.close()
	f.close()
	f = open(path +'InputTable/tab3_for_vary_eta_old.csv', 'r')
	g = open(path + 'InputTable/tab3_for_vary_eta.csv', 'r')
	c = 0
	for x, y in zip(f, g):
		c+=1
		if c == 1: continue
		yy = float(y.split(',')[-1])
		xx = float(x.split(',')[-1])
		yyy = y.split('.')[-1]

		if (abs(xx-yy)> 0.0001) or (len(yyy)>4): 
			print(' in: ' +  x)
			print('out: ' + y)

	g.close()
	f.close()




if __name__ == '__main__':
	main()










