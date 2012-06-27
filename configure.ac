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
        ,[NRINCL="-I$withval/include"]
        ,[AC_MSG_ERROR([Need to specify the NR path])])
AC_SUBST(NRINCL)


# Get the CosmoLib include path
AC_ARG_WITH(cosmolib,AS_HELP_STRING([--with-cosmolib]
        ,[location of CosmoLib])
        ,[COSMOLIBINCL="-I$withval/include"]
        ,[AC_MSG_ERROR([Need to specify the CosmoLib path])])
AC_SUBST(COSMOLIBINCL)

#Check is MOKA fits files handling is enabled
AC_ARG_ENABLE(moka,AS_HELP_STRING([--enable-moka]
        ,[Enable MOKA fits file handling])
        ,[okM=yes]
        ,[okM=no])
if test "$okM" = "yes"; then
   AC_DEFINE(WITH_MOKA,1,[Enable MOKA fits file handling])
fi

# Get the CFITS include path                                                                                      
AC_ARG_WITH(cfits,AS_HELP_STRING([--with-cfits]
        ,[location of cfits])
        ,[ok=yes]
        ,[ok=no])
if test "$ok" = "no"; then
   if test "$okM" = "no"; then
      CFITSINCL=""
   else
      AC_MSG_ERROR([Need to specify the cfits path])
   fi
else
   CFITSINCL="-I$withval/include"
fi
AC_SUBST(CFITSINCL)

# Get the CCFITS include path                                                                                       
AC_ARG_WITH(CCfits,AS_HELP_STRING([--with-CCfits]
        ,[location of CCfits])
        ,[ok=yes]
        ,[ok=no])
if test "$ok" = "no"; then
   if test "$okM" = "yes"; then
      AC_MSG_ERROR([Need to specify the CCfits path])
   else
      CCFITSINCL=""
   fi
else
   CCFITSINCL="-I$withval/include"
fi
AC_SUBST(CCFITSINCL)


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
