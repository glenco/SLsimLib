
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

#ifndef pi
#define pi  M_PI
#endif

#ifndef Grav
#define Grav  4.7788e-20  // G/c^2 in Mpc
#endif

#ifndef hplanck
#define hplanck  6.626068e-27  // in erg*sec
#endif

#ifndef inv_hplanck
#define inv_hplanck 1.50919e26 // in 1/(erg*sec)
#endif

#ifndef error_message
#define error_message
#define ERROR_MESSAGE() std::cout << "ERROR: file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef line_message
#define line_message
#define PRINT_LINE() std::cout << "file: " << __FILE__ << " line: " << __LINE__ << std::endl;
#endif

#ifndef bool_declare
#define bool_declare
typedef enum {FALSE, TRUE, MAYBE} Boo;
#endif

#ifndef lensquant_declare
#define lensquant_declare
typedef enum {dt,alpha,alpha1,alpha2,kappa,gamma1,gamma2,gamma3,invmag} LensingVariable;
#endif

#ifndef PosType_declare
#define PosType_declare
typedef double PosType;
#endif

#ifndef IndexType_declare
#define IndexType_declare
typedef unsigned long IndexType;
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
