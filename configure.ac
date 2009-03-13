AC_INIT("OpenMx", "0.1")
AC_PREREQ([2.61])
AC_CONFIG_AUX_DIR("build-aux")

# check the host type
AC_CANONICAL_HOST


# check for gcc compiler
AC_PROG_CC
AC_MSG_CHECKING(whether gcc is installed)
if test "$GCC" != yes; then
AC_MSG_RESULT(no)
AC_MSG_ERROR(gcc is not the default C compiler)
else
AC_MSG_RESULT(yes)
fi

# check for gcc version number
AX_GCC_VERSION

# check target architecture


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
