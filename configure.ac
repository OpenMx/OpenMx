AC_INIT("OpenMx", "0.1")
AC_PREREQ([2.61])
AC_CONFIG_AUX_DIR("inst/tools")

npsol_library='inst/npsol'

: ${R_HOME=`R RHOME`} 
if test -z "${R_HOME}"; then 
echo "could not determine R_HOME" 
exit 1 
fi 

# check for gcc compiler
AC_PROG_CC
AC_MSG_CHECKING(whether gcc is installed)
if test "$GCC" != yes; then
AC_MSG_RESULT(no)
AC_MSG_ERROR(gcc is not the default C compiler)
compiler_specific_flags=""
else
AC_MSG_RESULT(yes)
compiler_specific_flags="-Wall"
fi

# check the host type
AC_CANONICAL_HOST
case "${host_os}" in
linux-gnu*) 
    omx_platform=linux 
    arch_specific_linker_flags=""
    arch_specific_compiler_flags=""
    ;;
darwin*) 
    omx_platform=osx 
    arch_specific_linker_flags="-mmacosx-version-min=10.4"
    arch_specific_compiler_flags="-mmacosx-version-min=10.4"
    ;;
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
AC_SUBST(compiler_specific_flags)
AC_SUBST(arch_specific_linker_flags)
AC_SUBST(arch_specific_compiler_flags)
AC_CONFIG_FILES([src/Makevars])

AC_CONFIG_COMMANDS([src/omxSymbolTable.h],[${R_HOME}/bin/R --slave --vanilla < inst/tools/genSymbolTableHeader.R  > src/omxSymbolTable.h])
AC_CONFIG_COMMANDS([src/omxSymbolTable.c],[${R_HOME}/bin/R --slave --vanilla < inst/tools/genSymbolTableSource.R  > src/omxSymbolTable.c])
AC_OUTPUT
