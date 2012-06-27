#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([SLsimLib], [1.0], [margie.petkova@gmail.com])
AM_INIT_AUTOMAKE([1.10 no-define foreign])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


# Get the NR include path
AC_ARG_WITH(nr,AS_HELP_STRING([--with-nr]
        ,[location of nr])
        ,[NR_DIR="$withval"]
        ,[NR_DIR=""])
AC_SUBST(NR_DIR)


# Get the MOKAfits include path
AC_ARG_WITH(mokafits,AS_HELP_STRING([--with-mokafits]
        ,[location of MOKAfits])
        ,[MOKAFITS_DIR="$withval"]
        ,[MOKAFITS_DIR=""])
AC_SUBST(MOKAFITS_DIR)


# Get the CosmoLib include path
AC_ARG_WITH(cosmolib,AS_HELP_STRING([--with-cosmolib]
        ,[location of CosmoLib])
        ,[COSMOLIB_DIR="$withval"]
        ,[COSMOLIB_DIR=""])
AC_SUBST(COSMOLIB_DIR)

# Get the CFITS include path                                                                                      
AC_ARG_WITH(cfits,AS_HELP_STRING([--with-cfits]
        ,[location of cfits])
        ,[CFITS="$withval"]
        ,[CFITS="/usr/local/cfitsio"])
AC_SUBST(CFITS)

# Get the CCFITS include path                                                                                       
AC_ARG_WITH(CCfits,AS_HELP_STRING([--with-CCfits]
        ,[location of CCfits])
        ,[CCFITS="$withval"]
        ,[CCFITS="/usr/local/CCfits"])
AC_SUBST(CCFITS)


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
