source("_common.R")

library("care")
data(efron2004)
dim(efron2004$x) # 442 10
colnames(efron2004$x)

x_mat <- as.matrix(efron2004$x)   
class(x_mat) <- "matrix"        

X <- x_mat[,-2]         # remove "sex" (a binary variable)
y <- efron2004$y[,1]

fit_data <- list(
  N=nrow(X),
  y = y,
  P = ncol(X), 
  X = apply(X, 2, rescale_to_range1)
)

library(mgcv)

# Fit GAM with same structure as GP model: y = linear + smooth terms
# This matches the GP structure: y = beta0 + X*beta1 + f1(x1) + ... + f9(x9)

# Create data frame for GAM
gam_data <- data.frame(
  y = fit_data$y,
  fit_data$X  # X is already scaled from fit_data
)
# Make sure column names match the original predictor names
colnames(gam_data)[-1] <- colnames(efron2004$x)[-2]  # Remove "sex" column


### model 1: GAM with linear + smooth terms for each predictor
# s(xi, fx=FALSE) allows the smoother to determine its own complexity
gam_fit_time <- system.time({
  gam_fit <- gam(y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6 +
                     s(age) + s(bmi) + s(bp) + s(s1) + s(s2) + 
                     s(s3) + s(s4) + s(s5) + s(s6),
                 data = gam_data,
                 method = "REML")
})

# Summary of GAM fit
summary(gam_fit)

# Display time taken
gam_fit_time

# Plot GAM smooth functions
par(mar = c(3, 3, 2, 1), oma = c(0, 0, 0, 0))
plot(gam_fit, pages = 1, pch = 16, cex = 0.5, shade = TRUE, shade.col = "lightblue", bg = "red", col = "blue4", ylim=c(-2,2),lwd=2, seWithMean=TRUE)

## ===================== correct plot
# Set up a multi-panel plot
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(1, 1, 2, 1))

predictor_names <- c("age", "bmi", "bp", "s1", "s2", "s3", "s4", "s5", "s6")

for (pred in predictor_names) {
  # Plot original scatter points
  plot(gam_data[[pred]], gam_data$y, 
       pch = 16, cex = 0.5, col = "blue4",
       main = pred,
       xlab = pred, 
       ylab = "y")
  
  # Generate sequence for prediction
  x_seq <- seq(min(gam_data[[pred]]), max(gam_data[[pred]]), length.out = 100)
  
  # Create newdata for prediction: set current predictor to x_seq, others to their means
  newdata <- data.frame(matrix(rep(colMeans(gam_data[, -1]), each = 100), nrow = 100, byrow = FALSE))
  colnames(newdata) <- colnames(gam_data)[-1]
  newdata[[pred]] <- x_seq
  
  # Predict fitted values and SE
  pred_vals <- predict(gam_fit, newdata = newdata, se.fit = TRUE)
  
  # Add fitted line and confidence bands
  lines(x_seq, pred_vals$fit, col = "red", lwd = 2)
  lines(x_seq, pred_vals$fit + 1.96 * pred_vals$se.fit, col = "red", lty = 2)
  lines(x_seq, pred_vals$fit - 1.96 * pred_vals$se.fit, col = "red", lty = 2)
}

# Add overall title and reset plot layout
mtext("GAM Fitted Trends with Original Data Points", outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1))



### 选项2：强制线性+非线性分离
gam_fit_time <- system.time({
  gam_fit <- gam(y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6 +
                     s(age, fx=FALSE, bs="tp", by=I(age-mean(age))) + 
                     s(bmi, fx=FALSE, bs="tp", by=I(bmi-mean(bmi))) + 
                     s(bp, fx=FALSE, bs="tp", by=I(bp-mean(bp))) + 
                     s(s1, fx=FALSE, bs="tp", by=I(s1-mean(s1))) + 
                     s(s2, fx=FALSE, bs="tp", by=I(s2-mean(s2))) + 
                     s(s3, fx=FALSE, bs="tp", by=I(s3-mean(s3))) + 
                     s(s4, fx=FALSE, bs="tp", by=I(s4-mean(s4))) + 
                     s(s5, fx=FALSE, bs="tp", by=I(s5-mean(s5))) + 
                     s(s6, fx=FALSE, bs="tp", by=I(s6-mean(s6))),
                 data = gam_data,
                 method = "REML")
})

summary(gam_fit)

# Plot GAM smooth functions for separated linear + nonlinear model
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), oma = c(0, 0, 0, 0))

predictor_names <- c("age", "bmi", "bp", "s1", "s2", "s3", "s4", "s5", "s6")

for (i in 1:length(predictor_names)) {
  pred_name <- predictor_names[i]
  
  # 绘制原始数据点
  plot(gam_data[[pred_name]], gam_data$y, 
       pch = 16, cex = 0.5, col = "blue4",
       main = pred_name,
       xlab = pred_name, 
       ylab = "y")
  
  # 添加GAM拟合的平滑曲线
  # 注意：这里需要使用predict来获取完整的拟合值（线性+平滑）
  x_seq <- seq(min(gam_data[[pred_name]]), max(gam_data[[pred_name]]), length.out = 100)
  
  # 创建预测数据，其他变量设为均值
  newdata <- gam_data[1:100, ]  # 复制结构
  for (j in 1:length(predictor_names)) {
    if (j == i) {
      newdata[[predictor_names[j]]] <- x_seq
    } else {
      newdata[[predictor_names[j]]] <- mean(gam_data[[predictor_names[j]]])
    }
  }
  
  # 预测并绘制
  pred_vals <- predict(gam_fit, newdata = newdata, se.fit = TRUE)
  
  # 绘制拟合线和置信带
  lines(x_seq, pred_vals$fit, col = "red", lwd = 2)
  lines(x_seq, pred_vals$fit + 1.96 * pred_vals$se.fit, col = "red", lty = 2)
  lines(x_seq, pred_vals$fit - 1.96 * pred_vals$se.fit, col = "red", lty = 2)
}

### 选项2的简化版本：中心化平滑项
gam_fit_time <- system.time({
  # 先中心化数据
  gam_data_centered <- gam_data
  predictor_cols <- c("age", "bmi", "bp", "s1", "s2", "s3", "s4", "s5", "s6")
  
  # 为平滑项创建中心化变量
  for(col in predictor_cols) {
    gam_data_centered[[paste0(col, "_c")]] <- gam_data[[col]] - mean(gam_data[[col]])
  }
  
  gam_fit <- gam(y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6 +
                     s(age_c, bs="tp") + s(bmi_c, bs="tp") + s(bp_c, bs="tp") + 
                     s(s1_c, bs="tp") + s(s2_c, bs="tp") + s(s3_c, bs="tp") + 
                     s(s4_c, bs="tp") + s(s5_c, bs="tp") + s(s6_c, bs="tp"),
                 data = gam_data_centered,
                 method = "REML")
})




##### tow step model: linear first, then smooth on residuals
# Step 1: Fit a linear model
lm_fit <- lm(y ~ age + bmi + bp + s1 + s2 + s3 + s4 + s5 + s6, data = gam_data)

# Step 2: Compute residuals from the linear model
gam_data$residuals <- residuals(lm_fit)

# Step 3: Fit a GAM on the residuals
gam_fit_residuals <- gam(
  residuals ~ s(age, bs = "tp") + s(bmi, bs = "tp") + s(bp, bs = "tp") +
    s(s1, bs = "tp") + s(s2, bs = "tp") + s(s3, bs = "tp") +
    s(s4, bs = "tp") + s(s5, bs = "tp") + s(s6, bs = "tp"),
  data = gam_data,
  method = "REML"
)

# Step 4: Plot everything
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1), oma = c(1, 1, 2, 1))

predictor_names <- c("age", "bmi", "bp", "s1", "s2", "s3", "s4", "s5", "s6")

for (pred in predictor_names) {
  # Plot original scatter points
  plot(gam_data[[pred]], gam_data$y, 
       pch = 16, cex = 0.5, col = "blue4",
       main = pred,
       xlab = pred, 
       ylab = "y")
  
  # Generate sequence for prediction
  x_seq <- seq(min(gam_data[[pred]]), max(gam_data[[pred]]), length.out = 100)
  
  # Create newdata for prediction
  newdata <- data.frame(matrix(rep(colMeans(gam_data[, -1]), each = 100), nrow = 100, byrow = FALSE))
  colnames(newdata) <- colnames(gam_data)[-1]
  newdata[[pred]] <- x_seq
  
  # Predict linear trend
  linear_trend <- predict(lm_fit, newdata = newdata)
  
  # Predict smooth trend (from GAM on residuals)
  smooth_trend <- predict(gam_fit_residuals, newdata = newdata, type = "response")
  
  # Combine linear and smooth trends
  full_trend <- linear_trend + smooth_trend
  
  # Add linear trend (dashed green line)
  lines(x_seq, linear_trend, col = "darkgreen", lty = 2, lwd = 1.5)
  
  # Add full trend (solid red line)
  lines(x_seq, full_trend, col = "red", lwd = 2)
  
  # Add legend to the first plot
  if (pred == predictor_names[1]) {
    legend("topright", 
           legend = c("Full (Linear + Smooth)", "Linear only"),
           col = c("red", "darkgreen"),
           lty = c(1, 2),
           lwd = c(2, 1.5),
           cex = 0.7)
  }
}

# Add overall title and reset plot layout
mtext("Two-Step Model: Linear + Smooth Effects", outer = TRUE, cex = 1.2)
par(mfrow = c(1, 1))

# Summary of the linear model
summary(lm_fit)

# Summary of the GAM on residuals
summary(gam_fit_residuals)
