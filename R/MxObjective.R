setClass("MxObjective", representation("VIRTUAL"))

setConstructorS3("MxObjective", function() {

  extend(Object(), "MxObjective");

})

setMethodS3("createMxJobClosureR", "MxObjective", function(objective, job, ...) {});
setMethodS3("createMxJobClosureC", "MxObjective", function(objective, job, ...) {});