// Cut & paste from microbenchmark package (BSD license)

#include <inttypes.h>
typedef uint64_t nanotime_t;

#if defined(WIN32)
#include "nanotimer_windows.h"
#elif defined(u__MACH__) || defined(u__APPLE__)
#include "nanotimer_macosx.h"
#elif defined(linux) || defined(u__linux) || defined(u__FreeBSD__) || defined(u__OpenBSD__)
#include "nanotimer_gettime.h"
#elif defined(sun) || defined(u__sun) || defined(u_AIX)
#include "nanotimer_rtposix.h"
#else /* Unsupported OS */
nanotime_t get_nanotime(void) { return 0; }
#endif
