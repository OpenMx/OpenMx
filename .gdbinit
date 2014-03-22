dir src
set breakpoint pending on
b Rf_error
# _omxRaiseError is not necessarily an error; optimizer can test non-pd covariance matrices during search
b _omxRaiseError
catch throw
