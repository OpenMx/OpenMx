#
#   Copyright 2007-2016 The OpenMx Project
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


omxBrownie <- function(quantity=1, walnuts=TRUE){
	if (walnuts==FALSE)stop("Walnuts are required for brownies. Please correct the 'walnuts' argument to either TRUE (default) or 'allergic'.")
	amt <- c(4, 2, 8, 1.25, 1, 0.75, 0.5, 1.5)
	amt <- round(amt*quantity, 3)
	unit <- c("whole", "cup", "oz", "cup",
		"tsp", "cup", "tsp", "cup")
	ing <- c("eggs", "granulated sugar",
			"butter (room temperature)",
			"cocoa powder", "vanilla extract",
			"all-purpose flour", "salt",
			"chopped walnuts")
	brown <- list(
		matrix(c(amt, unit, ing), 
			ncol=3,
			dimnames=list(NULL, c("Qty.", "Unit", "Ingredient"))
			),
		matrix(c(quantity, quantity, 1, 1, 1, 1,
			"13 x 9 inch baking dish", 
			"24 inch piece of standard-width parchment paper",
			"Large metal or glass bowl",
			"Hand Mixer",
			"Spoon, Spatula or Other Similar Implement",
			"Toothpick",
			"Cooling Rack",
			"Standard Oven"), 
			ncol=2,
			dimnames=list(NULL, c("Qty.", "Object"))
			),
		matrix(c(
		"Preheat oven to 350 degrees Fahrenheit.", 
		"Lightly butter baking pan(s), then insert parchment sling.",
		"Beat eggs using mixer until homogonous. Beat in sugar until 'ribbon stage' is reached.",
		"Mix in other ingredients, stirring until just combined.",
		"Bake at 350 for 25 to 30 minutes. Brownies are done when toothpick inserted in brownies comes out clean.",
		"Place on rack to cool. When nearly at room temperature, brownies may be removed by picking up parchment sling for easy cutting."), ncol=1)
		)
	names(brown) <- c("Ingredients", "Equipment", "Procedure")
	if (walnuts=="allergic") {
		brown$Ingredients <- brown$Ingredients[1:7,]
	}
	return(brown)
}

### future work for OmxBrownie objects
### this doesn't work yet

# setClass(Class = "OmxBrownie",
# 	representation = representation(
# 		ingredients = "matrix", 
# 		equipment = "matrix", 
# 		procedure = "matrix",
# 		display = "character", "VIRTUAL"))
		
# setGeneric("imxCreateBrownie", 
# 	function(ingredients, equipment, procedure, ...) {
# 		return(standardGeneric("imxCreateBrownie"))
# })

# setMethod("imxCreateBrownie", "OmxBrownie",
# 	function(ingredients, equipment, procedure, ...) {
# 		.Object <- populateMatrixSlot(.Object, "labels", labels, nrow, ncol)
# 		.Object <- populateMatrixSlot(.Object, "values", values, nrow, ncol)
# 		.Object <- populateMatrixSlot(.Object, "free", free, nrow, ncol)
# 		.Object <- populateMatrixSlot(.Object, "lbound", lbound, nrow, ncol)
# 		.Object <- populateMatrixSlot(.Object, "ubound", ubound, nrow, ncol)
# 		.Object@name <- name
# 		return(.Object)
# 	}
# )	



