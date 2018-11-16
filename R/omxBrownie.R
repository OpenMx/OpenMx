#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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
	if ((is.character(walnuts) && walnuts == 'allergic') ||
		    (is.logical(walnuts) && walnuts == TRUE)) {
		if (is.character(walnuts)) walnuts <- FALSE
	} else {
		stop("Walnuts are required for brownies. Please correct the 'walnuts' argument to either TRUE (default) or 'allergic'.")
	}
	if (walnuts) {
		amt <- c(.5, .5, 4/3,   1, 1/4, 2/3, 3, 16, 180, 1, 1)
	} else {
		amt <- c(.5, .5, 2/3, 2/3, 1/4,   0, 3, 16, 180, 1, 1)
	}
	amt <- round(amt*quantity, 3)
	unit <- c("cup", "cup", "cup", "cup", "tsp",
		'cup', "cup", "dates", "grams", 'Tbsp', "tsp")
	ing <- c("whole wheat flour",
		"black bean flour",
		"sucanat",
		"Dutched cocoa powder",
		"salt",
		"chopped walnuts",
		"water",
		"deglet nour dates",
		"coconut butter (gently warmed to 90F)",
		"blackstrap molasses",
		"vanilla extract")
	brown <- list(
		data.frame("Qty."=amt, "Unit"=unit, "Ingredient"=ing),
		matrix(c(quantity, quantity, 1, 1, 1, 1, quantity, 1,
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
		"Heat oven to 350 degrees Fahrenheit.", 
		"Mix dry ingredients with a fork.",
		"Puree water and dates in a high speed blender.",
		"Combine all ingredients and mix thoroughly.",
		"Bake at 350 for 25 to 30 minutes. Brownies are done when toothpick inserted in brownies comes out clean.",
		"Place on rack to cool. When nearly at room temperature, brownies may be removed by picking up parchment sling for easy cutting."), ncol=1)
		)
	names(brown) <- c("Ingredients", "Equipment", "Procedure")
	if (!walnuts) {
		brown$Ingredients <- brown$Ingredients[-match("chopped walnuts", ing),]
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



