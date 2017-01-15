#include <windows.h>

nanotime_t get_nanotime(void) {
    LARGE_INTEGER time_var, frequency;
    QueryPerformanceCounter(&time_var);
    QueryPerformanceFrequency(&frequency);

    /* Convert to nanoseconds */
    return 1.0e9 * time_var.QuadPart / frequency.QuadPart;
}
