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


