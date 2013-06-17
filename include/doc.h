/** \file
 * doc.h
 *
 * \brief This is a file to define doxygen groups for documentation purposes
 *
 *  Created on: Oct 11, 2011
 *      Author: bmetcalf
 */

/** \mainpage
   <pre>
 GLAMER is a very flexible computer code for doing simulations of gravitational lensing 
 and modeling of gravitational lenses.
 
 This page describes  the steps to doing a simple simulation and where in the
 documentation to find further information.
 
   <em> GETTING STARTED: </em>
 
 GLAMER is in the form of a C++ library that can be linked into your code.  
 For instructions on installing and linking the library see 
 http://metcalf1.bo.astro.it/wiki/projects/glamer/GLAMER.html
 
   <em> READ IN PARAMETERS: </em>
 
 Within your main code the first thing to do is to read in the input parameter file.
 Your parameter file should contains all the parameters required to run the simulation.
 A sample parameter file is provided in the repository.   A particular simulation will
 not require all the parameters to be specified.  The parameter file should be read by 
 constructing an InputParams object.  See the ImputParams class for options.  If a parameter 
 is required, but not present in the parameter file the program should throw an exception 
 and notify you of which parameter needs to be included.
 
 The list of currently used input parameters is below.  They are also in the sample_paramfile
 that comes with the distribution.

		####### General ##########
		outputfile         output    # will contain the image magnifications

		deflection_off		0		# switches deflection off, default is 0
		
		read_redshift_planes   0 		# 0 - no, 1 - yes , reads in the redshifts of the lensing planes from a file
		redshift_planes_file    Z.txt   # an asci file with redshifts of the planes, excluding the source redshift

		####### Main halos type ##########
		main_halo_on	   1    # 0: no main halo, 1 or else - a main halo
		main_DM_halo_type       1    # DM internal profile type: 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
		main_galaxy_halo_type   1	# if set, a galaxy profile is chosen: 0 or none, 1 or NSIE
		
		###### Field halos ######
		field_off           0   # run without field halos, default is 0
		field_Nplanes       10  # number of field planes 
		
		###### Field halos type ##########

		field_internal_profile         1    # DM internal profile type: 0 or nolens, 1 or NFW, 2 or PseudoNFW, 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
		field_internal_profile_galaxy  1	# if set, a galaxy profile is chosen: 0 or none, 1 or NSIE
		field_galaxy_mass_fraction    0.1	# must be set if field_internal_profile_galaxy is set, mass fraction of the galaxy

		field_prof_internal_slope_pnfw		2 	# slope of the PseudoNFW profile, default is 2
		field_prof_internal_slope_pl		-1  # slope of the PowerLaw profile, default is -1

		###### Field halos from a mass function #######
		field_mass_func_type	   1		# mass function type, 0: Press-Schechter, 1: Sheth-Tormen, 2: PowerLaw
		field_fov		   1.0e4	        # field of view in square arcseconds
		field_buffer	   1.0		# in physical Mpc
		field_min_mass	   1.0e9			# min mass of the halos in the light cone in solar masses

		####### Field halos from an input file ##############
		#field_input_simulation_file ../MillenniumData/DMhalos.txt  # if set, the light cone is read from an input file

		####### AnaNSIE lens halo model ############
		sigma              250       # velocity dispersion in km/s
		core               0.0e-5    # core radius in Mpc
		axis_ratio         0.8       # axis ration
		pos_angle          0	     # inclination angle in degrees

		NDistortionModes   0	     # number of ellipsoid distortion modes
		beta_perturb       1.0
		kappa_peturb       0.03
		gamma_perturb      0.03
		monopole_perturb   0.0
		quadrapole_perturb 0.005
		hexopole_perturb   0.005
		octopole_perturb   0.01

		NdensitySubstruct   0.0e6    # number density of substructure
		beta_sub           -1.0
		alpha_sub          -1.9
		R_submax           0.5e-3
		sub_mass_max           1.0e9
		sub_mass_min           1.0e6
		sub_type           1

		Nstars             0	    # number of stars to be implanted
		fstars             0.50	    # stellar mass fraction
		stars_mass         0.5	    # star mass in solar masses

		######## Type of source SB model
		SourceSBType	   0	 # 0: Uniform, 1: Gaussian, 2: BLR_Disk, 3: BLR_Sph1, 4: BLR_Sph2

		###### Gaussian source model ###########
		gauss_r2	   0.0

		###### BLR source model ################
		BHmass             1.0e9        # BH mass of the quasar
		gamma              -0.5
		inclin             35.0
		opening_ang        10
		r_in               5.0e-9
		r_out              5.0e-6
		nuo                6.17284e14
		source_sigma       0

		####### General #############
		z_lens             0.42		# lens redshift
		z_source           3.62		# source redshift
 

   <em> CONSTRUCT MODEL: </em>
 
   <em> CONSTRUCT GRID: </em>
 
   <em> SHOOT RAYS: </em>
 
   <em> OUTPUT: </em>
 
   </pre>
 */

/** \defgroup HighLevel High Level Routines
 *
 * \brief These are routines that can be used without having to muck around in the internals.
 */


/** \defgroup ImageFinding Image Finding
 * \ingroup HighLevel
 *
 * \brief Routines for finding and characterizing images.
 *
 */

/** \defgroup ChangeLens Lens Manipulation
 * \ingroup HighLevel
 *
 * \brief Routines for adjusting the lens.
 *
 */

/** \defgroup FitLens Lens Fitting
 * \ingroup HighLevel
 *
 * \brief Routines for fitting a lens to observations.
 *
 */

/** \defgroup Image Image Processing
 * \ingroup HighLevel
 *
 * \brief These routines are used for putting images into pixels, simulating PSF smoothing and noise,
 * fitting profiles etc.
 *
 */


/** \defgroup MidLevel Mid Level Routines
 *
 * \brief These are routines that might be used by a user.
 */

/** \defgroup ConstructorL2 Constructors and Destructors
 * \ingroup MidLevel
 *
 * \brief Routines for initializing and freeing the data types.  These have not yet been moved to C++
 * form.  They a left over from C.
 *
 */

/** \defgroup ImageFindingL2 Image Finding
 * \ingroup MidLevel
 *
 * \brief Routines for finding and characterizing images.
 *
 */

/** \defgroup DeflectionL2 Deflection Solver
 * \ingroup MidLevel
 *
 * \brief Routines for calculating the deflection of a ray.
 *
 */

/** \defgroup FitLensL2 Lens Fitting
 * \ingroup MidLevel
 *
 * \brief Routines for fitting a lens to observations.
 *
 */

/** \defgroup ImageL2 Image Processing
 * \ingroup MidLevel
 *
 * \brief These routines are used for putting images into pixels, simulating PSF smoothing and noise,
 * fitting profiles etc.
 *
 */

/** \defgroup LowLevel Low Level Routines
 *
 * \brief These are routines that should never really be used by a user.
 */

/** \defgroup cosmolib Cosmology Library
 *
 * \brief Library for calculating all things cosmological (distances, power spectra, mass function, etc.)
 */

/** \defgroup Utill Utilities
 *
 * \brief Useful general purpose utilities.
 */

/** \defgroup function Functions
 *
 * \brief Global functions.  They are not all documented.
 *
 */

