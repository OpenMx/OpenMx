AC_DEFUN([AX_GCC_VERSION], [
  GCC_VERSION=""
  AX_GCC_OPTION([-dumpversion],[],[],[
    ax_gcc_version_option=yes
  ],[
    ax_gcc_version_option=no
  ])
  AS_IF([test "x$GCC" = "xyes"],[
    AS_IF([test "x$ax_gcc_version_option" != "no"],[
      AC_CACHE_CHECK([gcc version],[ax_cv_gcc_version],[
        ax_cv_gcc_version="`$CC -dumpversion`"
        AS_IF([test "x$ax_cv_gcc_version" = "x"],[
          ax_cv_gcc_version=""
        ])
      ])
      GCC_VERSION=$ax_cv_gcc_version
    ])
  ])
  AC_SUBST([GCC_VERSION])
])

AC_DEFUN([AX_GCC_OPTION], [
  AC_REQUIRE([AC_PROG_CC])

  AC_MSG_CHECKING([if gcc accepts $1 option])

  AS_IF([ test "x$GCC" = "xyes" ],[
    AS_IF([ test -z "$3" ],[
      ax_gcc_option_test="int main()
{
        return 0;
}"
    ],[
      ax_gcc_option_test="$3"
    ])

    # Dump the test program to file
    cat <<EOF > conftest.c
$ax_gcc_option_test
EOF

    # Dump back the file to the log, useful for debugging purposes
    AC_TRY_COMMAND(cat conftest.c 1>&AS_MESSAGE_LOG_FD)

    AS_IF([ AC_TRY_COMMAND($CC $2 $1 -c conftest.c 1>&AS_MESSAGE_LOG_FD) ],[
                AC_MSG_RESULT([yes])
        $4
    ],[
                AC_MSG_RESULT([no])
        $5
    ])
  ],[
    AC_MSG_RESULT([no gcc available])
  ])
])

AC_DEFUN([AX_C_CHECK_FLAG],[
  AC_PREREQ([2.61])
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AC_PROG_SED])

  flag=`echo "$1" | $SED 'y% .=/+-(){}<>:*,%_______________%'`

  AC_CACHE_CHECK([whether the C compiler accepts the $1 flag],
    [ax_cv_c_check_flag_$flag],[

    AC_LANG_PUSH([C])

    save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $1"
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([$2],[$3])
    ],[
      eval "ax_cv_c_check_flag_$flag=yes"
    ],[
      eval "ax_cv_c_check_flag_$flag=no"
    ])

    CFLAGS="$save_CFLAGS"

    AC_LANG_POP

  ])

  AS_IF([eval "test \"`echo '$ax_cv_c_check_flag_'$flag`\" = yes"],[
    :
    $4
  ],[
    :
    $5
  ])
])


