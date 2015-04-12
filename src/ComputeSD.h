#include <valarray>
#include <math.h>
#include "omxState.h"
#include "omxFitFunction.h"
#include "omxExportBackendState.h"
#include "Compute.h"
#include "matrix.h"

static void SD_grad(GradientOptimizerContext &, double);
void omxSD(GradientOptimizerContext &, int);

