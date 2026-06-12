library(OpenMx)
library(testthat)

# 1. Simulate raw data
set.seed(42)
n <- 2000
x1 <- rnorm(n)
x2 <- 0.5 * x1 + rnorm(n)
x3 <- 0.3 * x2 + rnorm(n)
x4 <- 0.7 * x3 + rnorm(n)
x5 <- 0.2 * x4 + rnorm(n)
data <- data.frame(x1, x2, x3, x4, x5)
manifests <- colnames(data)

# 2. Build RAM model template
make_model <- function() {
    mxModel("CI_benchmark_vector",
        type = "RAM",
        manifestVars = manifests,
        mxData(data, type = "raw"),
        mxPath(from = "x1", to = "x2", arrows = 1, free = TRUE, values = 0.5, labels = "beta_12"),
        mxPath(from = "x2", to = "x3", arrows = 1, free = TRUE, values = 0.3, labels = "beta_23"),
        mxPath(from = "x3", to = "x4", arrows = 1, free = TRUE, values = 0.7, labels = "beta_34"),
        mxPath(from = "x4", to = "x5", arrows = 1, free = TRUE, values = 0.2, labels = "beta_45"),
        mxPath(from = manifests, arrows = 2, free = TRUE, values = 1, labels = paste0("var_", manifests)),
        mxPath(from = "one", to = manifests, arrows = 1, free = TRUE, values = 0, labels = paste0("mean_", manifests)),
        mxFitFunctionML(vector = TRUE),
        mxCI(c("beta_12", "beta_23", "beta_34", "beta_45", "var_x1", "var_x2", "var_x3", "var_x4", "var_x5"))
    )
}

model_init <- make_model()
fit_init <- mxTryHard(model_init)

model_seq <- mxOption(fit_init, "Number of Threads", 1)
fit_seq <- mxRun(model_seq, intervals = TRUE)
ci_seq <- fit_seq$output$confidenceIntervals

model_par11 <- mxOption(fit_init, "Number of Threads", 11)

for (i in 1:10) {
    cat("Iteration", i, "running with 11 threads...\n")
    time_par11 <- system.time({
        fit_par11 <- mxRun(model_par11, intervals = TRUE)
    })
    ci_par11 <- fit_par11$output$confidenceIntervals
    
    # Check for NA
    num_nas <- sum(is.na(ci_par11))
    cat("  Time:", time_par11["elapsed"], "seconds. NAs found:", num_nas, "\n")
    
    # Check if identical
    tryCatch({
        expect_equal(ci_seq, ci_par11, tolerance = 1e-4)
        cat("  Success: Values are identical to baseline!\n")
    }, error = function(e) {
        cat("  FAILED: Mismatch found!\n")
        print(ci_par11)
    })
}
