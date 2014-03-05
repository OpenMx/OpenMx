dir src
set breakpoint pending on
b Rf_error
b _omxRaiseError
catch throw
