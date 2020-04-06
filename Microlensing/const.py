import numpy as np
from MLG.const import * 
__all__ = ['const_Einsteinradius']
const_Einsteinradius  = 180/ np.pi * 3.6e6 * np.sqrt(4 * G_Newton * M_sun / (c_light ** 2 * pcinm * 1000))   #(mas/Msun)**0.5
