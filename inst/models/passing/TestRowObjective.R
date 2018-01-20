#
#   Copyright 2007-2018 The OpenMx Project
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


#-----------------------------------------
# Author: Michael Hunter
# Filename: TestRowObjective.R
# Purpose: Have a very basic test of the 
#  mxFitFunctionRow function
#-----------------------------------------

#-----------------------------------------
# Revision History
# 	Date: Wed Jun 9 15:57:04 EDT 2010 - Created Test
# 	Date: Fri Sep 10 14:33:13 EDT 2010
# 	Date: Mon Oct 4 16:29:10 EDT 2010 - Added omxCheckCloseEnough
#										Added to models/failing
#										Added OpenMx License
#
#-----------------------------------------


library('OpenMx')

set.seed(14) # Make repeatable

xdat <- data.frame(a=rnorm(10), b=1:10) # Make data set

# Make model
rmod <- mxModel(
		name = 'EZRowObjectiveTest',
		mxData(observed=xdat, type='raw'),
		mxAlgebra(sum(filteredDataRow), name='rowAlgebra'),
		mxAlgebra(sum(rowResults), name='reduceAlgebra'),
		mxFitFunctionRow('rowAlgebra', 'reduceAlgebra', dimnames = c('a', 'b'))
)
#rmodns <- mxOption(rmod, 'No Sort Data', c('a', 'b'))
#rmodnsRun <- mxRun(rmodns)
rmodRun <- mxRun(rmod)
# summary(rmodRun) #So far I haven't got the summary to work
# rmod$algebras


#-----------------------------------------
# Check that the row algebra works
#  Row algebras are evaluated row-wise
#  Results of the row-wise evaluation are stored in the row results

hardAdd <- mxEval(rowResults, rmodRun)
# 'G' should be the same as
ezAdd <- as.matrix(xdat$a + xdat$b, ncol=1)

omxCheckCloseEnough(hardAdd, ezAdd, 10^(-6))


#-----------------------------------------
# Check that the reduce algebra works
#  Reduce algebras combine the elements of the row results

hardTotal <- mxEval(reduceAlgebra, rmodRun)
# 'R' should be the same as
ezTotal <- sum( xdat$a + xdat$b )

omxCheckCloseEnough(hardTotal, ezTotal, 10^(-6))

