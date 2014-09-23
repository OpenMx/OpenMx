# R --no-save -f IntroSEM-UnivariateRaw.R --args 1
require(OpenMx)
library(httr)

args <- commandArgs(trailingOnly = TRUE)

host <- '127.0.0.1'
port <- 1337
server <- paste(host,port,sep=':')
apiurl <- paste("http://", server, "/api", sep="")

# need something cryptographically secure?
name <- paste0("agent", sample.int(1e7, 1))

data(multiData1)
parts <- cut(1:nrow(multiData1), 4)  # chop into 4 partitions
mask <- parts == levels(parts)[ as.integer(args[[1]]) ]

message("Waiting for model to be published")

while (1) {
  r <- try(GET(paste0(apiurl, "/model")), silent = TRUE)
  if (!is(r, "try-error") &&
        any(unlist(content(r)$models) == "test")) break  # wait for our model to appear
  Sys.sleep(1)
}

message("Found model")

r <- GET(paste0(apiurl, "/model/test"))
uniRegModelRaw <- unserialize(charToRaw(content(r)$model))

uniRegModelRaw <- mxModel(uniRegModelRaw,
                          mxData(observed=multiData1[mask,], type="raw"),
                          mxComputeOnce('fitfunction', 'fit'))

parNames <- names(omxGetParameters(uniRegModelRaw, FALSE, NA))

evaluation <- -1

while (1) {
  r <- GET(paste0(apiurl, "/model/test/param"))
  cr <- content(r)
  
  if (cr$evaluation < 1 || cr$evaluation == evaluation) {
    # We already addressed these parameter vectors
    Sys.sleep(1)
    next
  }
  evaluation <- cr$evaluation
  
  par <- unlist(cr$at)
  uniRegModelRaw <- omxSetParameters(uniRegModelRaw, labels=parNames, values=par)
  uniRegModelRawOut <- mxRun(uniRegModelRaw, silent = TRUE)
  
  r <- POST(paste0(apiurl, "/model/test/fit"), encode="json",
            body=list(agent=name, evaluation=cr$evaluation, fit=uniRegModelRawOut$output$fit))
  Sys.sleep(1)
}
