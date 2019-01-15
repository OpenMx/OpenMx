# avoid user interaction
options(repos = c(CRAN = "https://cloud.r-project.org"))

ap = available.packages()
pk = c('ctsem', 'EasyMx', 'ifaTools', 'metaSEM', 'semtree', 'umx')
v = sapply(pk, function(x) ap[which(ap[,'Package'] == x), 'Version'],
            USE.NAMES=FALSE)
cat(mapply(function(p1, v1) { paste0(p1,'_',v1,'.tar.gz') },
           pk, v, USE.NAMES = FALSE))
