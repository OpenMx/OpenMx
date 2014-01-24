#
#   Copyright 2007-2012 The OpenMx Project
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


# Taken from the package "parallel" available in R >= 2.14.0
omxDetectCores <- function(all.tests = FALSE, logical = FALSE) {
	if("package:parallel" %in% search()) {
		return(detectCores(all.tests, logical))
	} else if(.Platform$OS.type == "windows") {        
            return(1L)
    } else {
		systems <-
			list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
			freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
			linux = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
			irix  = c("hinv | grep Processors | sed 's: .*::'",
						"hinv | grep '^Processor '| wc -l"),
			solaris = if(logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/usr/sbin/psrinfo -p")
		for (i in seq(systems)) {
        	if(all.tests ||
            	length(grep(paste("^", names(systems)[i], sep=''),
                               R.version$os))) {
                for (cmd in systems[i]) {
                	a <- gsub("^ +","", system(cmd, TRUE)[1])
                    if (length(grep("^[1-9]", a))) return(as.integer(a))
                }
			}
		}
        return(NA_integer_)
    }
}
