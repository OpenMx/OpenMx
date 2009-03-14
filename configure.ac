AC_INIT("OpenMx", "0.1")
AC_PREREQ([2.61])
AC_CONFIG_AUX_DIR("build-aux")

npsol_library='inst/npsol'

# check for gcc compiler
AC_PROG_CC
AC_MSG_CHECKING(whether gcc is installed)
if test "$GCC" != yes; then
AC_MSG_RESULT(no)
AC_MSG_ERROR(gcc is not the default C compiler)
else
AC_MSG_RESULT(yes)
fi

# check the host type
AC_CANONICAL_HOST
case "${host_os}" in
linux-gnu*) omx_platform=linux ;;
darwin*) omx_platform=osx ;;
*) AC_MSG_ERROR(host os is not recognized: ${host_os}) ;;
esac

case "${omx_platform}" in
linux)
   case "${host_cpu}" in
   *86) omx_arch=x86 ;;
   *86_64) omx_arch=x86_64 ;;
   *) AC_MSG_ERROR(host arch is not recognized: ${host_cpu}) ;;
   esac

   # check for gcc version number
   AX_GCC_VERSION
   AC_PROG_AWK
   major_release=`echo "$GCC_VERSION" | $AWK -F '.' '{print $1 "." $2}'`
   npsol_library="${npsol_library}/linux/${omx_arch}/gcc${major_release}"
;;
osx)
   npsol_library="${npsol_library}/osx"
;;
*)
   AC_MSG_ERROR(openmx platform is not recognized: ${omx_platform})
esac

npsol_file="${npsol_library}/libnpsol.a"
AC_CHECK_FILE([${npsol_file}],,[AC_MSG_ERROR([npsol library not found])])

npsol_dir="../${npsol_library}"
AC_SUBST(npsol_dir)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
