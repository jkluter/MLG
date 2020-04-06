import warnings
import os
warnings.filterwarnings('ignore')

#Version 
Version = '1'
Subversion = '0'
Date = '01.04.2020'
print('MLG Version %s.%s %s' % (Version, Subversion, Date))

#defult parameter
path = __file__[:-11]
inputpath =  path + 'InputTable/'
default_table = inputpath + 'events_within_gaia.fits'
imagepath = path + 'Images/'
paperpath =  path[:-4]+ 'Simulation_paper/' 
message_folder = os.path.expanduser('~')+ '/Dropbox/Data/'

from . import Simulation
from . import StellarMotion
from . import Modeling
from . import Math
from . import Microlensing
from . import Analysis
from . import name
from . import const
from . import image
from . import utils


