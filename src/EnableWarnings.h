#ifdef ENABLEWARNINGS_H_
#error "Only include EnableWarnings.h once from a *.cpp file"
#else
#define ENABLEWARNINGS_H_
#endif

#ifdef EXTRA_GCC_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#pragma GCC diagnostic warning "-Wmisleading-indentation"
#pragma GCC diagnostic warning "-Wextra"
#pragma GCC diagnostic warning "-Wunused-function"
#pragma GCC diagnostic warning "-Wunused-local-typedefs"
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic warning "-Wignored-attributes"
#endif
