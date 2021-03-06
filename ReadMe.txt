MLG
MicroLensing with Gaia 


Requirements:

Python3.7.3
numpy 1.16.2
scipy 1.2.1
joblib 0.13.2
matplotlib 3.0.3
astropy 3.1.2
imageio 2.6.1


Usage:

In order to use this code copy the folder to any local directory.

The code can be used as either as package: 

 
>>> import MLG

or as executable file

$ python MLG -keywords

in the local directory

The simulation is controneld either via keywords (see below) or via a  GUI wich is started if no keywords are set. 

This GUI can be used to
1. Start the Simulation:   Press "Start"
2. Print an comand for starting the simulation: Press "Print command"
3. Starting the Analysis GUI: Press "Analysis"

For starting the simulation the following keywords can be used: 

list of events (-j):  

	Single events: (-j 1)
		Use default table "InputTable/events_within_gaia.fits"
		each Row consist of one Lens-Source pair
	Multiple events: 
		Use default table "InputTable/1.fits" ... "InputTable/22.fits"
		each Table consinst of one Lens and multiple Sources
		the following selection of sources are possible

		All sources:  (-j 2 / -j 4/ -j 6 / -j 8)   
			use all sources for the simulation 
		Sources with pm: (-j 3/ -j 4 / -j 7/ -j 8  )
			use only sources with a 5 paramerter solution in Gaia Dr2 
		Sources with sigma< 0.5: (-j 5/ -j 6 / -j 7 / -j 8 )   
				select all sources with a 5 parameter soluion in Gaia Dr2 and 
				an expectet precission in along-scan direction smaller than 0.5"

		Any combination of the three cases can be startet at once
		Be aware of long computaion times for "All sources" 

	Other: -j location_to_file
		use the given table with one microlensing event per line.
		The table must be readable by astropy.tabel.Table and 
		 must be contained the following colums:
		
		Gaia source id of the Lens: 'source_id'
		Gaia source id of the Source: 'ob_source_id'

		Position of the Lens in deg: 'ra' / 'dec'
		Position of the Source  in deg: 'ob_ra' / 'ob_dec'
		Error of the position of the Lens in deg: 'ra_error' / 'dec_error'
		Error of the position of the Source  in deg: 'ob_ra_error' / 'ob_dec_error'

		Proper motion of the Lens in mass/year: 'pmra' / 'pmdec'
		Proper motion of the source in mass/year: 'ob_pmra' / 'ob_pmdec'
		Error of the proper motion of the Lens in mass/year: 'pmra_error' / 'pmdec_error'
		Error of the proper motion of the source in mass/year: 'ob_pmra_error' / 'ob_pmdec_error'

		Parallax of the Lens in mas: 'parallax'
		Parallax of the Source in mas: 'ob_parallax'
		Error of the parallax of the Lens in mas: 'parallax_error'
		Error of the parallax of the Source in mas:'ob_parallax_error'

		G Magnitude of the Lens: 'phot_g_mean_mag'
		G Magnitude of the Source: 'ob_phot_g_mean_mag'
		Mass of the lens in solar mass: 'mass'
		Error of the ,ass in solar mass: 'mass_error'
		Epoch of the closest approach in julian year: 'tca'

		It might be needed to include observations of the events extracted from GOST 
		(https://gaia.esac.esa.int/gost/) in the "InputTable/Observations.csv" and 
		"InputTable/Observations_extended.csv" tables
File string: -s 
	Identifier for the result file.
	Multiple string can be used separated by a " " (The GUI also except ",")

Observations:
	Use of the extended 10 years mission: -e 
	Use of N external 2D observations: -obs N  (N = 2 if N is not given )

Monte-Carlo-Simulation
	N Picks from the error ellips: -ne N (default N = 500 )
	Vary imput parameters: -np N (default N = 100 )
		use "-v 0" if imput parameters should not be varied 
Parallel computing:
	Number of threads: -c N (default N = 6) 


Output: 
	one or more pkl files cotaining the fited parameters for each realisiation
	can be loaded by 
		% MLG.Simulation.MonteCarlo.loadPkl(File)
		Contains:
			 	'Header': List of all control imputs, and Code Version
			 	'Results': List of the fit parameters
			 	'Input_parameter': List of the input parameters
			 	'Results_ext_obs':  List of the fit parameters with external_obs if ext_obs > 0 
			 	'Chi2':	List of the reduced CHI2
			 	'Eventlist': (see Input)
			 	'Names': (see Input)






