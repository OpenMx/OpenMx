/* 
 * Solaris does not define gethrtime if _POSIX_C_SOURCE is defined because
 * gethrtime() is an extension to the POSIX.1.2001 standard. According to [1]
 * we need to define __EXTENSIONS__ before including sys/time.h in order to
 * force the declaration of non-standard functions.
 *
 * [1] http://www.oracle.com/technetwork/articles/servers-storage-dev/standardheaderfiles-453865.html
 */
#if defined(sun) || defined(__sun)
#define __EXTENSIONS__
#endif

# include <sys/time.h>

/* short an sweet! */
nanotime_t get_nanotime(void) {
    return gethrtime();
}
