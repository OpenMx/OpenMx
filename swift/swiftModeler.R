# generates all models 
# write out as mxobjects
# process all objects

connections <- 4
numcol      <- noquote(strsplit(allinputs," ")[[1]][2])
rowcol 	    <- as.integer(numcol)
size = rowcol*rowcol

# creates the .adat files 

### call a swift job with the array of possible 
### number of connections and total size
### swift will generate the model objects in parallel

system("swift -tc.file tc.data -sites.file sites_ncsa modelgen.swift -user=skenny -connections=4") 

### resulting model objects will be in modelObjects subdir
### call swift script for processing the model objects
### with simplecovariane

system("swift -tc.file tc.data -sites.file sites_ncsa modelproc.swift -user=skenny")

