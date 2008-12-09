#
# Objective Functions are for optimizing!
#
setClass(Class = "MxFIMLObjective",
	representation = representation(
		covariance = "MxSymmetricMatrix",
		means = "numeric"))
		

