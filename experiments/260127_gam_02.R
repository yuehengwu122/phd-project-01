# Load required libraries
library(gam)  # For Generalized Additive Models (GAM)
library(mgcv) # For modern GAM fitting with penalized splines

# -----------------------------------------------------------------------------
# Diabetes Dataset
# -----------------------------------------------------------------------------

# Load and preprocess the diabetes dataset
data(efron2004, package = "care")
x_mat <- as.matrix(efron2004$x)   # Convert predictors to matrix
class(x_mat) <- "matrix"         # Ensure class is "matrix"
X <- x_mat[,-2]                   # Remove "sex" (a binary variable)

# Create data frame for GAM analysis
gam_data <- data.frame(
  y = efron2004$y[,1],
  apply(X, 2, rescale_to_range1)
)
colnames(gam_data)[-1] <- colnames(efron2004$x)[-2]  # Match predictor names

# Fit GAM using the `gam` package
gam_diabetes <- gam::gam(
  y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6 +
      s(age) + s(bmi) + s(bp) + s(s1) + s(s2) + 
      s(s3) + s(s4) + s(s5) + s(s6),
  data = gam_data
)
summary(gam_diabetes)

# Fit GAM using the `mgcv` package with thin-plate splines
preds <- c("age","bmi","bp","s1","s2","s3","s4","s5","s6")
form_full <- as.formula(paste0(
  "y ~ ",
  paste(preds, collapse = " + "),
  " + ",
  paste0("s(", preds, ", bs='tp', m=1, k=10)", collapse = " + ")
))
mgcv_diabetes <- mgcv::gam(form_full, data = gam_data, method = "REML")
summary(mgcv_diabetes)

# -----------------------------------------------------------------------------
# Boston Housing Dataset
# -----------------------------------------------------------------------------

# Load and preprocess the Boston Housing dataset
data("Boston", package = "MASS")
gam_data <- data.frame(
  y = Boston$medv,  # Median value of owner-occupied homes
  as.matrix(Boston[,-c(4,9,14)])  # Exclude categorical predictors
)

# Fit GAM using the `gam` package
gam_boston <- gam::gam(
  y ~ crim + zn + indus + nox + rm + age + dis + tax + ptratio + black + lstat +
      s(crim) + s(zn) + s(indus) + s(nox) + s(rm) + s(age) + s(dis) + s(tax) + s(ptratio) + s(black) + s(lstat),
  data = gam_data
)
summary(gam_boston)
anova(gam_boston)  # ANOVA for parametric effects

# Fit GAM using the `mgcv` package with thin-plate splines
preds <- c("crim","zn","indus","nox","rm","age","dis","tax","ptratio","black","lstat")
form_full <- as.formula(paste0(
  "y ~ ",
  paste(preds, collapse = " + "),
  " + ",
  paste0("s(", preds, ", bs='tp', m=1, k=10)", collapse = " + ")
))
mgcv_boston <- gam(form_full, data = gam_data, method = "REML")
summary(mgcv_boston)

