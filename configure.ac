#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([cnvCaller],[1.4.0],[yannis.duffourd@u-bourgogne.fr])

# Store the auxiliary build tools (e.g., install-sh, config.sub, config.guess)
# in this dir (build-aux)
# AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_AUX_DIR([.])

# init automake 
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Safety checks in case user overwritten --srcdir
# AC_CONFIG_SRCDIR([cnvCompare.cpp utils.cpp main.cpp])

AC_CONFIG_MACRO_DIR([m4])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# add optimisation flag 
# AX_APPEND_FLAG([-O3])

# check support of ISO C++20
# AX_CHECK_COMPILE_FLAG(-std=c++2a, AX_APPEND_FLAG([-std=c++2a]), AC_MSG_ERROR([-std=c++2a not supported, you may have to set CXXFLAGS to enable C++2a support.]))
AX_CHECK_COMPILE_FLAG(-std=c++2a, AX_APPEND_FLAG([-std=c++2a]), AC_MSG_ERROR(["std=c++2a not supported, you may have to set CXXFLAGS to enable C++2a support."]))

# Checks for libraries.
# BOOST
AC_CHECK_LIB([boost_program_options], , AC_MSG_ERROR([Unable to find the Boost Program option library]))





# Checks for header files.
AC_CHECK_HEADERS([sys/time.h])
AC_CHECK_HEADERS([sys/types.h])
AC_CHECK_HEADERS([boost/program_options.hpp])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
