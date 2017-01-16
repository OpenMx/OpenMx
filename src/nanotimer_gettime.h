#include <time.h>

static const nanotime_t nanoseconds_in_second = 1000000000LL;

nanotime_t get_nanotime(void) {
    struct timespec time_var;

    /* Possible other values we could have used are CLOCK_MONOTONIC,
     * which is takes longer to retrieve and CLOCK_PROCESS_CPUTIME_ID
     * which, if I understand it correctly, would require the R
     * process to be bound to one core.
     */
    clock_gettime(CLOCK_MONOTONIC, &time_var);

    nanotime_t sec = time_var.tv_sec;
    nanotime_t nsec = time_var.tv_nsec;

    /* Combine both values to one nanoseconds value */
    return (nanoseconds_in_second * sec) + nsec;
}
