#include <mach/mach_time.h>

/* see http://developer.apple.com/library/mac/#qa/qa2004/qa1398.html */
nanotime_t get_nanotime(void) {
    uint64_t time;
    mach_timebase_info_data_t info;
    
    time = mach_absolute_time();
    mach_timebase_info(&info);
    
    /* Convert to nanoseconds */
    return time * (info.numer / info.denom);
}
