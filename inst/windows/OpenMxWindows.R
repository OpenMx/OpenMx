install.packages('snowfall')
winDialog("ok", "Now select the OpenMx zip file when prompted")
filename <- file.choose()
install.packages(filename, repos=NULL)