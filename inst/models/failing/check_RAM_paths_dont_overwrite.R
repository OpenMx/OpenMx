## This behavior is purposeful.  We made a design decision to not allow mxPath 
## to create simultaneous single headed and double headed arrows in the same 
## position in the A and S matrices.  The rare (really, are there any?) cases 
## when this makes sense are far outweighed by the times when people are either 
## (a) making a mistake or (b) trying to convert a single headed into a double 
## headed arrow or vice versa.  Of course, this behavior can be overridden by 
## directly manipulating an mxMatrix.  We did (and I do still) consider mxPath() 
## to be the training wheels for OpenMx.  The idea being we would like to make 
## it difficult to make mistakes.  I believe this philosophy is still correct.
##
## https://lists.virginia.edu/sympa/arc/openmx-developers/2016-03/msg00085.html

# "~/bin/OpenMx/inst/models/failing/check_RAM_paths_dont_overwrite.R"

# When adding paths in mxModel, I think everything the user requests should be added
# in line with what you say is what you get?

m1 <- mxModel("x_to_y_is_OK", type="RAM",
	manifestVars = c("X", "Y"),
	# both these paths should be added:
	mxPath("X", to = "Y", arrows = 1)
)
# X->Y is fine
omxCheckEquals(m1$A$free["Y","X"], TRUE)

m1 <- mxModel("x_with_y_overwrites_x_to_y", type="RAM",
	manifestVars = c("X","Y"),
	# both these paths should be added:
	mxPath("X", to = "Y", arrows = 1), 
	mxPath("X", to = "Y", arrows = 2)
)
# Oops: the 2-headed path overwrote the 1-headed path
omxCheckEquals(m1$A$free["Y","X"], TRUE)

# Last-path-in is OK
omxCheckEquals(m1$S$free["Y","X"], TRUE)


