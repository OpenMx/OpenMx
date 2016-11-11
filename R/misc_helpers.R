#' omx_make_fake_data
#'
#' This function takes as argument an existing dataset, which 
#' must be either a matrix or a data frame. Each column of the 
#' dataset must consist either of numeric variables or ordered 
#' factors. When one or more ordered factors are included, 
#' then a heterogeneous correlation matrix is computed using 
#' John Fox's polycor package. Pairwise complete observations 
#' are used for all covariances, and the exact pattern of 
#' missing data present in the input is placed in the output,
#' provided a new sample size is not requested. Warnings from
#' the polycor::hetcor function are suppressed.
#'
#' @param dataset The original dataset of which to make a simulacrum
#' @param digits = Round the data to the requested digits (default = 2)
#' @param n Number of rows to generate (NA = all rows in dataset)
#' @param use.names Whether to name the variables (default = TRUE)
#' @param use.levels = Whether to use existing levels (default = TRUE)
#' @param use.miss Whether to have data missting as in original (defaults to TRUE)
#' @param mvt.method = Passed to hetcor (default = "eigen")
#' @param het.ML = Passed to hetcor (default = FALSE)
#' @param het.suppress Passed to hetcor (default = TRUE)
#' @return - new dataframe
#' @family Data Functions
#' @export
#' @examples
#' fakeCars = omx_make_fake_data(mtcars)
omx_make_fake_data <- function(dataset, digits = 2, n = NA, use.names = TRUE, use.levels = TRUE, use.miss = TRUE, mvt.method = "eigen", het.ML = FALSE, het.suppress = TRUE){
  # requires mvtnorm & polycor
  # requires data frame or matrix
  if((is.data.frame(dataset)+is.matrix(dataset)) == 0){
    warning("Data must be a data frame or matrix")
  }
  # organization
  row <- dim(dataset)[1] # number of rows
  if(is.na(n))(n <- row) # sets unspecified sample size to num rows
  col <- dim(dataset)[2] # number of columns
  del <- is.na(dataset)  # records position of NAs in dataset
  if(n != row){
    select <- round(runif(n, 0.5, row + .49),0)
    del    <- del[select,]
  }
  num <- rep(NA, col)    # see what's not a factor
  ord <- rep(NA, col)    # see what's an ordered factor

  # which columns are numeric (the others are factors)?
  for (i in 1:col){
    num[i] <- is.numeric(dataset[,i])
    ord[i] <- is.ordered(dataset[,i])
  }

  # check for unordered factors
  location <- !(num|ord)
  unorder  <- sum(location)

  if(unorder>0)warning(
    paste("Unordered factor detected in variable(s):", 
      names(dataset)[location]
    )
  )

  # if everything is numeric, don't invoke polycor
  if(sum(!num) == 0){
    # generate data with rmvnorm
	# depends on mvtnorm::rmvnorm
    fake <- mvtnorm::rmvnorm(n, apply(dataset, 2, mean, na.rm = TRUE),
		cov(dataset, use = "pairwise.complete.obs"), mvt.method)

    # round the data to the requested digits
    fake <- round(fake, digits)

    # insert the missing data, if so requested
    if(use.miss == TRUE)(fake[del] <- NA)

    # give the variables names, if so requested
    if(use.names == TRUE)(names(fake) <- names(dataset))

    # return the new data
    return(fake)
  }

  # if there are factors, we start here

  # find the variable means (constrain to zero for factors)
  mixedMeans <- rep(0, col)
  mixedMeans[num] <- apply(dataset[, num], 2, mean, na.rm = TRUE)

  # estimate a heterogeneous correlation matrix
  if (het.suppress == TRUE){
	  suppressWarnings(het <- polycor::hetcor(dataset, ML = het.ML))
  } else {
	  het <- polycor::hetcor(dataset, ML = het.ML)	
  }
  mixedCov <- het$correlations

  # make a diagonal matrix of standard deviations to turn the 
  # correlation matrix into a covariance matrix
  stand <- matrix(0, col, col)
  diag(stand) <- rep(1, col)
  diag(stand)[num] <- apply(dataset[,num], 2, sd, na.rm = TRUE)
  # pre and post multiply hetero cor matrix by diagonal sd matrix
  mixedCov <- stand %*% mixedCov %*% stand

  # generate the data
  fake <- as.data.frame(mvtnorm::rmvnorm(row, mixedMeans, mixedCov, mvt.method))

  # insert the missing data, if so requested
  if(use.miss == TRUE)(fake[del] <- NA)

  # turn the required continuous variables into factors
  for (i in (1:col)[!num]){
    # the original data for this column
    old <- dataset[,i]
   
    # the new data for this column, omiting NAs
    new <- fake[!is.na(fake[,i]),i]

    # what are the levels of the original factor?
    lev <- levels(old)

    # establish cutpoints in new variable from cdf of old factor
    cut <- cumsum(table(old))/(sum(!is.na(old)))

    # put continuous variable into a matrix, repeating value across columns
    wide <- matrix(new, length(new), length(lev))

    # put the cutpoints in a matrix, repeating the cut point values across rows
    crit <- matrix(quantile(new, cut), length(new), length(lev), byrow = TRUE)

    # for each value (row of the wide matrix), 
    # how many cutpoints is the value greater than?
    # number of cutpoints surpassed=category
    fake[!is.na(fake[,i]),i] <- apply(wide>crit, 1, sum)

    # make it a factor
    fake[,i] <- factor(fake[,i], ordered = TRUE)

    # give the new factor the same levels as the old variable
    if(length(levels(fake[,i]))!=length(lev))message(
      paste("Fewer categories in simulated variable", 
      names(fake)[i], "than in input variable", names(dataset)[i]))
	if(use.levels == TRUE & (length(levels(fake[,i])) == length(lev))){
		levels(fake[,i]) <- lev
  	} else {
		levels(fake[,i]) <- 1:length(lev)
	}
  }

  # Round the data to the requested digits
  fake[,num] <- round(fake[,num], digits)

  # Give the variables names, if so requested
  if(use.names == TRUE)(names(fake) <- names(dataset))
  
  # Return the new data
  return(fake)
}