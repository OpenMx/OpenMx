library(OpenMx)
library(httr)

host <- '127.0.0.1'
port <- 1337
server <- paste(host,port,sep=':')
apiurl <- paste("http://", server, "/api", sep="")

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

r <- POST(paste0(apiurl, "/model/test"), body=list(model=blob), encode="json")

#r <- POST(paste0(apiurl, "/model/rabit"), body=list(model=blob), encode="json")

par <- omxGetParameters(uniRegModelRaw)

uniRegModelRaw <- mxModel(
  mxMatrix(values=par, labels=names(par), nrow=length(par), ncol=1, free=TRUE),
  mxComputeGradientDescent(),
  mxFitFunctionR(fitfun = function(model, state) {
  par <- omxGetParameters(model)
  r <- PUT(paste0(apiurl, "/model/test/param"), body=list(param=par), encode="json")
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

expectSE <-c(0.053849, 0.071873, 0.104204, 0.047674, 0.078154)

expectMin <- 3151.492

omxCheckCloseEnough(expectVal, uniRegModelRawOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
                    as.vector(uniRegModelRawOut$output$standardError), 0.001)

omxCheckCloseEnough(expectMin, uniRegModelRawOut$output$minimum, 0.001)


