#
#   Copyright 2007-2009 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

beginTravisFold <- function(key, text) {
  cat(paste0("travis_fold:start:",key,"\033[33;1m",text,"\033[0m"),fill=T)
}

endTravisFold <- function(key) {
  cat(paste0("\ntravis_fold:end:",key,"\r"),fill=T)
}

args <- commandArgs(trailingOnly = TRUE)
if (any(args == 'gctorture')) {
	gctorture(TRUE)
	cat("*** GCTORTURE ENABLED ***\n");
}

if (all(file.exists("inst/models", "models") == c(TRUE,FALSE))) {
	file.symlink("inst/models", "models")
}

rexec <- Sys.getenv("REXEC")
if (rexec == "") rexec <- "R"
if (system2(rexec, args='--vanilla', stdout=FALSE, stderr=FALSE,
	    input='library(OpenMx); .Call(OpenMx:::.EigenDebuggingEnabled); q(status=3)') == 3) {
	cat("**", fill=TRUE)
	cat("** WARNING: Eigen conformability checking is DISABLED **", fill=TRUE)
	cat("**", fill=TRUE)
}

library(OpenMx)

options('mxPrintUnitTests' = FALSE)

directories <- c('models/passing')

if (any(args == 'nightly')) {
	directories <- c(directories, 'demo', 'models/nightly')
}
if (any(args == 'failing')) {
	directories <- c('models/failing')
}

null <- tryCatch(suppressWarnings(file('/dev/null', 'w')),  
	error = function(e) { file('nul', 'w') } )


if (any(args == 'gctorture')) {
	files <- c('models/passing/AlgebraComputePassing.R',
		   'models/passing/TestRowObjective.R')
} else {
	files <- list.files(directories, pattern = '^.+[.]R$',
			    full.names = TRUE, recursive = TRUE)
}

if (any(args == 'lisrel')) {
	files <- grep("LISREL", files, value=TRUE, ignore.case=TRUE)
}

errors <- list()
warnRec <- list()
runtimes <- numeric()

errorRecover <- function(script, opt, index) {
	mxSetDefaultOptions()
	mxOption(NULL, "Default optimizer", opt)
	sink(type = 'output')
	cat(paste(opt, index, "of",
		length(files), script, "...\n"))
	sink(null, type = 'output')
	start <- Sys.time()

	tryCatch.W.E <- function(expr) {
    W <- list()
	  w.handler <- function(w) { # warning handler
	    W[[1 + length(W)]] <<- w
	    invokeRestart("muffleWarning")
	  }
	  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler), warning = W)
	}
	
	got <- tryCatch.W.E(source(script, chdir = TRUE))
  
	stop.tm <- Sys.time()
	timeDifference <- stop.tm - start
	runtimes[[paste(opt,script,sep=":")]] <<- as.double(timeDifference, units = "secs")
  
	err <- got$value
	if (is(err, "error") && err$message != 'SKIP') {
	  errors[[opt]][[script]] <<- err$message
	  sink(type = 'output')
	  cat("*** ERROR from", script, '***\n')
	  print(err$message)
	  sink(null, type = 'output')
	}
	warnRec[[opt]][[script]] <<- got$warning
	
	rm(envir=globalenv(), 
		list=setdiff(ls(envir=globalenv()), 
			     c('warnRec', 'errors', 'errorRecover', 'opt', 'null',
			       'files', 'directories', 'runtimes', 'beginTravisFold', 'endTravisFold')))
}

optimizers <- strsplit(Sys.getenv("IMX_OPT_ENGINE"), "\\s")[[1]]

if (length(optimizers) == 0) {
	optimizers <- c('CSOLNP')
	if (!any(args == 'gctorture')) {
		optimizers <- c(optimizers, 'SLSQP')
		if (imxHasNPSOL()) {
			optimizers <- c(optimizers, 'NPSOL')
		}
	}
}

for (opt in optimizers) {
	beginTravisFold(opt, paste("TEST", opt))
	errors[[opt]] <- list()
	warnRec[[opt]] <- list()
	if (length(files) > 0) {
		for (i in 1:length(files)) {
			errorRecover(files[[i]], opt, i)
		}
	}
	sink(type = 'output')
	endTravisFold(opt)
}	

close(null)

totalErrors <- sum(sapply(errors, length))
cat("Number of errors:", totalErrors, '\n')
if (totalErrors > 0) {
	for (opt in names(errors)) {
		oerr <- errors[[opt]]
		fileName <- names(oerr)
		if (length(oerr)) for (i in 1:length(oerr)) {
			cat("Error", opt, fileName[[i]], '***\n')
			print(oerr[[i]])
			cat('\n')
		}
	}
} else {
	for (opt in names(warnRec)) {
		owarn <- warnRec[[opt]]
		count <- sapply(owarn, length)
    if (any(count > 0)) {
      cat("**", opt ,"total warnings =", sum(count), "\n", fill=TRUE)
      wscript <- names(owarn[count > 0])
      for (ws1 in wscript) {
        cat("*", opt, "warnings from", ws1, "\n", fill=TRUE)
        wlist <- owarn[[ws1]]
        for (w1 in wlist) {
          cat(paste("  ", w1), fill=TRUE)
        }
      }
    }
	}
}

write.csv(as.data.frame(runtimes), "runtimes.csv")

cat("Finished testing models.\n")
quit(status=totalErrors)
