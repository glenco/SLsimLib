Installing GLAMER                                                     {#install}
=================

Installing *GLAMER* is a mostly straightforward process, thanks to the *CMake*
build system which handles configuration options, dependency resolution, and
integration into projects.


The build system
----------------

The *GLAMER* libraries are built using the [CMake] build system. *CMake* is a
meta build system which generates the files required by your build environment
of choice, such as

-   Unix Makefiles
-   Xcode
-   Eclipse
-   ...


### Getting CMake

The installation of *CMake* on all major platforms is straightforward. There are
binary installation packages available on the *CMake* website.

For Linux, *CMake* can be found in most package managers natively.

For Mac OS X, *CMake* is available from both [Homebrew], [MacPorts] or [conda].


### Using CMake

*CMake* is a command-line tool that takes input files (called `CMakeLists.txt`)
and processes them into files for your build system of choice (via a generator),
optionally taking into account project-specific options. The normal workflow is
the following:

1.  Get source code containing `CMakeLists.txt` files for cmake to process.
2.  Create a new subdirectory which will contain the generated files. This is,
    in the language of cmake, an "out-of-source" build and the preferred mode
    of operation, since it does not contaminate the source directory (and cmake
    generates quite a number of files). The subdirectory is usually called
    "build".
3.  Call cmake with appropriate arguments to generate files for your build
    system (or IDE).
4.  Build using your usual method of choice.


Dependancies
------------

Before starting to configure and build *GLAMER*, it is advisable to resolve all
necessary dependencies, since the cmake configuration will automatically detect
settings based on the libraries it can find.

All the dependancies should be available through the package managment system [conda].  It is recommended that you use this to install and update the dependencies.

The *GLAMER* library requires some libraries by default and others are optional. 


Library   | Required by   || Description
----------|---------------||-----------------------------------------------------
[CFITSIO] | `ENABLE_FITS` |required| Library for reading and writing files in the FITS format.
[CCfits]¹ | `ENABLE_FITS` |required| C++ wrapper for *CFITSIO*.
[FFTW]    | `ENABLE_FFTW` |required| Library for Fast Fourier Transform algorithms.
[GSL]     | `ENABLE_GSL`  |optional| Library for scientific computation
[HDF5 C++]    | `ENABLE_HDF5` |optional| I/O in hdf5 format


Building the library
--------------------

After extracting the library, you will be find a directory structure that
contains folders for the individual *GLAMER* sub-libraries, together with an
additional `examples` folder containing a sample project, and a `cmake` folder
containing necessary files for *CMake*.

In order to build *GLAMER* using cmake, the following steps are necessary:

1.  Create a `build` directory for building the *GLAMER* libraries.
    
        $ mkdir build
        $ cd build
    
2.  Run `cmake` with arguments. You will probably want to set the generator and
    maybe some of the *GLAMER*-specific options. Here we generate Unix Makefiles
    and enable *FITS* and *FFTW* functions.
    
        $ cmake .. -G "Unix Makefiles" [-DENABLE_GSL=ON]
    
3.  At this point, you are presented with a Makefile or project inside the
    `build` directory which you can use to build the library.
    
        $ make

The last step of course depends on the type of generator you choose (using the
`-G` flag). For IDEs like Xcode or Eclipse, you have to open the project that
*CMake* generated and build the library from within the IDE.


### Rebuilding

A very useful feature of cmake is that it caches all the configuration options
you set in previous runs. Thus it is not necessary to specify options you want
to keep unchanged on subsequent runs. Continuing the above example, you could
turn off *FFTW* support using

    $ cmake .. -DENABLE_FFTW=OFF

while keeping the generator and `ENABLE_FITS` options untouched.

It is necessary to regenerate the build files when sources are added or removed.
Luckily, your build system should automatically check whether `cmake` needs to
be rerun for you before building.


Build options
-------------

Options can be set from the command-line using

    cmake .. -G <generator> -D<option>=<value>

For a list of available generators, use `cmake -h`.

The following options can be passed to `cmake` for building *GLAMER*:

Option        | Values      | Default | Description
--------------|-------------|---------|-----------------------------------------
`USE_CXX11`   | `ON`, `OFF` | `ON`    | Use the C++11 standard if supported by the compiler.
`N_THREADS`   | 1, 2, ...   | 1       | Set to the maximum number of threads to be used.
`ENABLE_FITS` | `ON`, `OFF` | `ON`    | Enable functions that need FITS support.
`ENABLE_FFTW` | `ON`, `OFF` | `ON`    | Enable functions that need FFTW support.
`ENABLE_GSL`  | `ON`, `OFF` | `OFF`   | Enable functions that need GSL support.
`_OPENMP`     | `ON`, `OFF` | `OFF`   | Enables openMP multi-threading.

More detailed descriptions of the individual options can be found below.

The list of currently cached option values can be shown using

    $ cmake .. -L


### Option `ENABLE_GSL`

This enables some halo model calculations in CosmoLib.


### Option `_OPENMP`

This is used in only one place and should be considered obsolete.


Using GLAMER
------------

Once the *GLAMER* libraries are built inside your build tree, they can be linked
with your executables as usual. However, you might want to try and manage your
project with *CMake* instead.

Having built the *GLAMER* libraries with CMake makes it easy to integrate them
into projects. Most of the process is automated by *CMake*. It is important to
have built *GLAMER* before setting up a project, since *CMake* then remembers
where to find the libraries.


### Creating new projects

The easiest way to start a project is to copy the sample project from the
`examples/project` folder or clone one of the example projects from the [Glenco Github page]("https://github.com/glenco"), specifically [Example1]("https://github.com/glenco/Example1"), [ExampleImage]("https://github.com/glenco/ExampleImage") and [ParticleExample]("https://github.com/glenco/ParticleExample").  The new project should be in its own directory outside of the *GLAMER* directory tree.


In the enclosed `CMakeLists.txt` file, all instances
of the "sample" name have to be replaced with the actual name of the project.
Sources and headers can be added to the appropriate lists.


### Existing cmake projects

For projects already using *CMake*, integrating *GLAMER* is as simple as calling

    find_package(GLAMER NO_MODULE REQUIRED)

before creating any targets. This will load the *GLAMER* configuration, import
the library targets and set variables containing the include directories and
libraries necessary for *GLAMER*, including all dependencies. These can then be
set using the usual

    include_directories(${GLAMER_INCLUDE_DIRS})

and

    target_link_libraries(<target> ${GLAMER_LIBRARIES})

functions.


### Building projects

Projects are built using *CMake* in the same way the libraries are built. After
creating a build tree, cmake will configure and generate a project which can
subsequently be used with the build tool of choice.

    $ cd myproject
    $ ls
    CMakeLists.txt  main.cpp
    $ mkdir build
    $ cmake .. -G Xcode
    $ open myproject.xcodeproj
or

    $ cd myproject
    $ ls
    CMakeLists.txt  main.cpp
    $ mkdir build
    $ cmake .. 
    $ make
If you are not using Xcode.


**Note:**
Building the project will not build the GLAMER libraries automatically!


[cmake]: http://www.cmake.org "CMake"
[homebrew]: http://www.brew.sh "Homebrew — The missing package manager for OS X"
[macports]: https://www.macports.org "The MacPorts Project"
[cfitsio]: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html "CFITSIO"
[ccfits]: https://heasarc.gsfc.nasa.gov/fitsio/CCfits/ "CCfits"
[ccfitsdoc]: https://heasarc.gsfc.nasa.gov/fitsio/CCfits/html/index.html "CCfits documentation"
[fftw]: http://www.fftw.org "FFTW Home Page"
[gsl]: http://www.gnu.org/software/gsl/ "GNU Scientific Library"
[HDF5 C++]: https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html"
[conda]:https://docs.conda.io/projects/conda/en/latest/"