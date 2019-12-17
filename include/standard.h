
#ifndef _STANDARD_DECLARE_
#define _STANDARD_DECLARE_

//#define NDEBUG  // Un-commenting this line will remove all assert() statements in the executable.

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cfloat>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pthread.h>

#include <list>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <limits>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <valarray>

#include <nr.h>
#include <nrD.h>
#include <utilities.h>

#include <cosmo.h>
#include <halo.h>

#ifndef PI
#define PI  M_PI
#endif

#ifndef Grav
#define Grav  4.7788e-20  // G/c^2 in Mpc
#endif

#ifndef arcsecTOradians
#define arcsecTOradians  0.000004848136811  // convert arcesconds to radians
#endif

#ifndef degreesTOradians
#define degreesTOradians  0.01745329251994  // convert degrees to radians
#endif

#ifndef hplanck
#define hplanck  6.626068e-27  // in erg*sec
#endif

#ifndef inv_hplanck
#define inv_hplanck 1.50919e26 // in 1/(erg*sec)
#endif

#ifndef kmpersecTOmpcperday
#define kmpersecTOmpcperday 2.80003e-15
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef line_message
#define line_message
#define PRINT_LINE() std::cout << "file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef boo_declare
#define boo_declare
typedef enum {NO, YES, MAYBE} Boo;
#endif

#ifndef lensquant_declare
#define lensquant_declare
enum LensingVariable {DELAYT,ALPHA,ALPHA1,ALPHA2,KAPPA,GAMMA,GAMMA1,GAMMA2,GAMMA3,INVMAG,PHI} ;
#endif

#ifndef PosType_declare
#define PosType_declare
typedef double PosType;
#endif

#ifndef IndexType_declare
#define IndexType_declare
typedef size_t IndexType;
#endif

#ifndef swap_declare
#define swap_declare
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#endif

#ifndef treeNBdim
#define treeNBdim 2  // dimension of boxes in tree
#endif

#ifdef DOUBLE_PRECISION
typedef double KappaType;
#else
typedef float KappaType;
#endif

// unit test definitions
// GlamerTest overrides these as necessary
#ifndef GLAMER_TEST
#define GLAMER_TEST_USES(t)
#define GLAMER_TEST_FRIEND(t)
#endif

#endif


#ifndef MpcToSeconds
#define MpcToSeconds (3.085677e22 / 2.99792458e8)
// Computed as MpcToMeter(m) / SpeedOfLight(m.s^{-1})
#endif

#ifndef SecondToDays
#define SecondToDays (1. / 86400.)
#endif

#ifndef SecondToYears
#define SecondToYears (1. / 86400. / 365.25)
// Taking 365.25 days per year.
#endif






