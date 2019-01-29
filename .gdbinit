set print thread-events off
dir src
set breakpoint pending on
# Rf_errorcall is used for some out of memory conditions
b Rf_errorcall
# _omxRaiseError is not necessarily an error; optimizer can test non-pd covariance matrices during search
#b _omxRaiseError
# FitContext::recordIterationError is used for soft feasibility constraints
#b FitContext::recordIterationError
catch throw
# Hints:
# Use Rf_PrintValue(SEXP) to pretty print R data from gdb
# Use pda(mem, rows, cols) to dump a matrix of doubles in column major order
# Use pia(mem, rows, cols) to dump a matrix of ints in column major order
