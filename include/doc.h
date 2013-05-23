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
 
 The list of currently used input parameters are:
 
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

