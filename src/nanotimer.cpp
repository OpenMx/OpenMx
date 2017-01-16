// Cut & paste from microbenchmark package (BSD license)

#include <inttypes.h>
typedef uint64_t nanotime_t;

#if defined(WIN32)
#include "nanotimer_windows.h"
#elif defined(__MACH__) || defined(__APPLE__)
#include "nanotimer_macosx.h"
#elif defined(linux) || defined(__linux) || defined(__FreeBSD__) || defined(__OpenBSD__)
#include "nanotimer_gettime.h"
#elif defined(sun) || defined(__sun) || defined(_AIX)
#include "nanotimer_rtposix.h"
#else /* Unsupported OS */
nanotime_t get_nanotime(void) { return 0; }
#endif
