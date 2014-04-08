set print thread-events off
dir src
set breakpoint pending on
b Rf_error
# _omxRaiseError is not necessarily an error; optimizer can test non-pd covariance matrices during search
b _omxRaiseError
catch throw
# Hints:
# Use Rf_PrintValue(SEXP) to pretty print R data from gdb
# Use pda(mem, rows, cols) to dump a matrix of doubles in column major order
# Use pia(mem, rows, cols) to dump a matrix of ints in column major order
