# Generate some sample data
set.seed(123)
x <- 1:100
y <- c(rnorm(50, 2*x[1:50], sd = 10), rnorm(50, 3*x[51:100], sd = 10))
plot(x, y)

# Fit a piecewise linear regression with two segments
library(segmented)
seg.model <- segmented(lm(y ~ x), seg.Z = ~x)
seg.model
lines(x, predict(seg.model))


# Calculate the variance of the intersection point
var1 <- var(predict(seg.model)[x <= seg.model$psi])
var2 <- var(predict(seg.model)[x > seg.model$psi])
var.intersect <- var1 + var2

# Calculate the predicted coefficient of variation
cv.predicted <- sqrt(var.intersect) / mean(x)

# Calculate the confidence interval for the intersection point
alpha <- 0.05  # significance level
df <- length(x) - 4  # degrees of freedom
t.val <- qt(1 - alpha/2, df)  # t-value for two-tailed test
se <- cv.predicted * mean(x) / t.val  # standard error
ci.low <- seg.model$psi - se
ci.high <- seg.model$psi + se
cat("The 95% confidence interval for the intersection point is",
    "[", ci.low, ",", ci.high, "].")


# for constructing the 95% CI based on original paper:
# estimate the variance in the y-values by using the residuals?
# var(resid - predicted)

# Chat GPT says

# Now, if we assume that the errors in the regression lines are normally distributed, then we can use the properties of a bivariate normal distribution to estimate the variance of the x-location of the intersection point x0. Specifically, we can use the fact that the correlation between the predicted y-values for each regression line is equal to the cosine of the angle between the regression lines.
#
# Using some trigonometry, we can show that the cosine of the angle between the regression lines is equal to (b1 - b2) / sqrt(1 + (b1 * b2)), where b1 and b2 are the slopes of the two regression lines.

# # Calculate the slope of each segment of the piecewise linear regression model
# b1 <- coef(lm_left)[1]
# b2 <- coef(lm_right)[2]
#
# # Estimate the variance in the x-location of the intersection point
# var.x0 <- (var1 * var2) / ((b1 - b2)^2 * (1 + (b1 * b2)))

# add
