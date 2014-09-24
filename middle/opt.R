library(OpenMx)
library(httr)
library(jsonlite)

if (packageVersion("httr") < package_version("0.5.0.9000")) {
  stop("A newer version of httr is required. You may need to install from https://github.com/hadley/httr")
}

host <- '127.0.0.1'
port <- 1337
server <- paste(host,port,sep=':')
apiurl <- paste("http://", server, "/api", sep="")

partitions <- 4

#r <- GET(apiurl)
#content(r)

manifests <- c("x1", "y")

uniRegModelRaw <- mxModel("FIML Univariate Regression of y on x1",
                          type="RAM",
                          manifestVars=manifests,
                          mxPath(from="x1", to="y", arrows=1, 
                                 free=TRUE, values=.2, labels="b1"),
                          mxPath(from=manifests, 
                                 arrows=2, free=TRUE, values=.8, 
                                 labels=c("VarX1", "VarE")),
                          mxPath(from="one", to=manifests, 
                                 arrows=1, free=TRUE, values=.1, 
                                 labels=c("MeanX1", "MeanY")))

blob <- rawToChar(serialize(connection=NULL, uniRegModelRaw, ascii=TRUE))

r <- POST(paste0(apiurl, "/model/test"),
          content_type_json(),
          body=toJSON(list(model=blob, minRows=partitions), digits=8))

#r <- POST(paste0(apiurl, "/model/rabit"), body=list(model=blob), encode="json")

par <- omxGetParameters(uniRegModelRaw)

uniRegModelRaw <- mxModel(
  mxMatrix(values=par, labels=names(par), nrow=length(par), ncol=1, free=TRUE),
  mxComputeGradientDescent(),
  mxFitFunctionR(fitfun = function(model, state) {
    par <- omxGetParameters(model)
    r <- PUT(paste0(apiurl, "/model/test/param"),
             body=toJSON(list(param=par), digits=8),
             content_type_json())
    #content(r)
    
    fit <- NULL
    while (1) {
      Sys.sleep(1)
      r <- GET(paste0(apiurl, "/model/test/fit"))
      if (!is.null(content(r)$nan)) {
        fit <- NA
        break
      }
      fit <- content(r)$fit
      if (length(fit)) break
    }
    print(fit)
    return(fit)
  }))

uniRegModelRawOut <- mxRun(uniRegModelRaw, suppressWarnings=TRUE)

#---------------------
# check values: uniRegModelRawOut

expectVal <- c(0.669179, 1.13643, 1.647629, 0.984894, 3.189368)

expectMin <- 3151.492 / partitions

omxCheckCloseEnough(expectVal, uniRegModelRawOut$output$estimate, 0.01)

omxCheckCloseEnough(expectMin, uniRegModelRawOut$output$minimum, 0.001)
