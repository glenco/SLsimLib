Documentation                                                        {#mainpage}
=============

*GLAMER* is a flexible computer code for doing simulations of gravitational
lensing and modeling of gravitational lenses.

This page describes the steps to doing a simple simulation and where in the
documentation to find further information. The installation and usage of
*GLAMER* is described on the [installation page](@ref install). If you would
like to obtain the code, please see the [appropriate section](@ref copying).

The code can be made available to you through GitHub.com upon request.  Requests should be sent to Ben Metcalf at robertbenton.metcalf@unibo.it.

Getting started
---------------

`GLAMER` is in the form of a C++ library that can be linked into your code.  
For instructions on installing and linking the library see [here](http://metcalf1.bo.astro.it/wiki/projects/glamer/GLAMER.html)


Read in Parameters
----------------------------

Original all objects where constructed from a parameter file through the InputParams class.  Many of these constructors still exist, but this method has been depricated.  It is recommended that the classes be constructed from thier other constructors now.



Original all objects where constructed from a parameter file through the InputParams class.  Many of these constructors still exist, but this method has been depricated.  It is recommended that the classes be constructed from thier other constructors now.

LensHalo classes
-----------------
To represent the mass that couses the lensing one or more LensHalo objects are constructed.  A full list of LensHalos can be found [here](file:///Users/bmetcalf/Glamer/Doc/html/classes.html).  Some examples are:



`LensHaloNFW` - a NFW halo

`LensHaloMassMap` - a mass sheet represented with a pixelized map

`LensHaloParticles` - a mass distribution represented by particles or point masses possibly from a simulation

`LensHaloRealNSIE` - a non-singular isothermal ellipical lens

`LensHaloPowerLaw` - a lens with a power-law mass profile



Source Class
------------------

Sources represent anything that produces a visible image infront of, within or behind the lens.  

Some of them are:

`SourceUniform` : Makes the source a uniform surface brightness circle

`SourceGaussian` : Source with a Gaussian profile

`SourceShapelet` : a source made from a shapelet decomposition 

`SourceSersic` : a source with a Sersic profile

`SourcePixelled` : a source made from a pixelized image of an object


Lens class
----------

The lens class contains the LensHalos and keeps the information about the light-cone.  

LensHalos are inserted into the Lens object with there redshift and angular position.  Once it is created the Lens owns the LensHalos. 

`Lens::ResetSourcePlane()` is used to change the redshift of the source plane.  It can be at lower redshift than some or all of the LensHalos.

Main halos are ones that are inserted indiviually. Field halos are those that are generated as a population throughout the light cone.  These are stored differently within the Lens.

Grid & GridMap
--------------

The `Grid` structure contains the image points and thier source points along with other information. Without further grid refinement the Grid can be used to make a map of deflection, convergence, shear or time-delay. Grid contains functions for outputing these. If no further grid refinement is desired for image or caustic finding, the `GridMap` structure can be used which requires signifcantly less memory overhead.

To change the source plane surface brightness use `Grid::RefreshSurfaceBrightnesses ()`.



Image finding and Grid refinement
---------------------------------

The mean functions used for image finding are

-   `ImageFinding::find_images_kist()` 
-   `ImageFinding::map_images_fixedgrid` and
-   `ImageFinding::map_images()`.
And for finding critical curves / caustics
-   `ImageFinding::find_crit()`

These take a `Grid` and a `Lens` object. The found images are then stored in an
array of `ImageInfo` structures.


Observational effects 
----------------------

The `Observation` class is used to simulate pdf, noise, etc.


Output
------

The `Grid` has some direct output functions for making maps that do not depend
on sources (kappa, gamma, etc.). The ImageInfo structure has some information
about the images found and contains a linked list of all the points within the
image. Each point contains information about the kappa, gamma, source position,
time-delay and in some cases the the surface brightness at that point.

