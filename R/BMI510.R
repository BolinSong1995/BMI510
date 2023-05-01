# bmi510.R
library(pwr)


#' Randomly sample elements or rows from a vector or data frame
#'
#' @param x An atomic vector, data frame, or tibble
#' @param n n samples will be sampled from X
#' @param replace default TRUE, this allows sampling with replacement
#'
#' @return vector or data frame
#' @export
rando = function(x, n = 1, replace = T) {
  if (is.data.frame(x) || is.matrix(x)) {
    return(x[sample(nrow(x), n, replace = replace), ])
  } else {
    return(sample(x, n, replace = replace))
  }
}


#' Check if values in a vector are equal to the minimum value
#'
#' This function accepts an atomic vector \code{x} and returns a logical vector with \code{TRUE} values where \code{x} equals its minimum value.
#'
#' @param x An atomic vector to check
#' @param na.rm Logical. Should missing values be removed? (default = TRUE)
#'
#' @return A logical vector with \code{TRUE} values where \code{x} equals its minimum value
#'
#' @examples
#' is_min(c(3, 5, 1, 3, 2))
#' is_min(c(3, 5, NA, 3, 2))
#'
#' @export
is_min = function(x, na.rm = T) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  x == min(x)
}

#' Check if values in a vector are equal to the max value
#'
#' This function accepts an atomic vector \code{x} and returns a logical vector with \code{TRUE} values where \code{x} equals its max value.
#'
#' @param x An atomic vector to check
#' @param na.rm Logical. Should missing values be removed? (default = TRUE)
#'
#' @return A logical vector with \code{TRUE} values where \code{x} equals its maximum value
#'
#' @examples
#' is_max(c(3, 5, 1, 3, 2))
#' is_max(c(3, 5, NA, 3, 2))
#'
#' @export
is_max = function(x, na.rm = T) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  x == max(x)
}


#' repeat matrix
#'
#' @param x A matrix or data frame to replicate
#' @param M The number of times to repeat the rows of \code{x} (default = 1)
#' @param N The number of times to repeat the columns of \code{x} (default = 1)
#'
#' @return A matrix with rows and/or columns replicated M and/or N times
#'
#' @examples
#' rep_mat(mtcars, 2, 1)
#' rep_mat(1:5, 3, 2)
#'
#' @export
rep_mat = function(x, M = 1, N = 1) {
  if(sum(class(X) %in% c(c("tbl_df","tbl","data.frame")))>0 | is.matrix(X)==T){
    temp = X[rep(1:nrow(X), M), ]
    temp = temp[, rep(1:ncol(temp), N)]
    return(temp)
  }else{
    stop("X is not data frame or matrix")
  }
}


#' This function takes a tibble \code{x} and returns a character vector containing the classes of each variable in the tibble.
#'
#' @param x A tibble to get the variable classes from
#'
#' @return A character vector containing the classes of each variable in the tibble
#'
#' @examples
#' classes(mtcars)
#'
#' @import tibble
#' @export
classes = function(x) {
  return(sapply(x, class))
}

#' Returns a tibble x in which the numeric variables have been scaled with scale
#'
#' @param x A tibble to scale
#' @param center Logical. Should the variables be centered? (default = TRUE)
#' @param scale Logical. Should the variables be scaled to have unit variance? (default = TRUE)
#'
#' @return A tibble with the numeric variables centered and/or scaled
#'
#' @examples
#' df_scale(mtcars)
#'
#' @importFrom tibble tibble
#' @importFrom stats scale
#' @export
df_scale = function(x, center = T, scale = T) {
  numeric_vars <- sapply(x, is.numeric)
  if (center) {
    x[numeric_vars] <- lapply(x[numeric_vars], scale, center = center, scale = FALSE)
  }
  if (scale) {
    x[numeric_vars] <- lapply(x[numeric_vars], scale, center = FALSE, scale = scale)
  }
  return(x)
}


#' This function calculates the log-likelihood of a sample \code{x} under the normal distribution, with mean \code{mean} and standard deviation \code{sd}.
#'
#' @param x A numeric vector of data
#' @param mean The mean of the normal distribution
#' @param sd The standard deviation of the normal distribution
#'
#' @return The log-likelihood of \code{x} under the normal distribution
#'
#' @examples
#' log_likelihood_norm(rnorm(100), 0, 1)
#'
#' @importFrom stats dnorm
#' @export
log_likelihood_norm = function(x, mean, sd) {
  return(sum(dnorm(x, mean, sd, log = T)))
}


#' This function calculates the log-likelihood of a sample \code{x} under the uniform distribution, with minimum value \code{min} and maximum value \code{max}.
#'
#' @param x A numeric vector of data
#' @param min The minimum value of the uniform distribution
#' @param max The maximum value of the uniform distribution
#'
#' @return The log-likelihood of \code{x} under the uniform distribution
#'
#' @examples
#' log_likelihood_unif(runif(100), 0, 1)
#'
#' @importFrom stats dunif
#' @export
log_likelihood_unif = function(x, min, max) {
  return(sum(dunif(x, min, max, log = T)))
}


#' This function calculates the log-likelihood of a sample \code{x} under the chi-squared distribution, with degrees of freedom \code{df}.
#'
#' @param x A numeric vector of data
#' @param df The degrees of freedom of the chi-squared distribution
#'
#' @return The log-likelihood of \code{x} under the chi-squared distribution
#'
#' @examples
#' log_likelihood_chisq(rchisq(100, 5), 5)
#'
#' @importFrom stats dchisq
#' @export
log_likelihood_chisq = function(x, df) {
  return(sum(dchisq(x, df, log = T)))
}


#' This function calculates the log-likelihood of a sample \code{x} under the F distribution, with degrees of freedom \code{df1} and \code{df2}.
#'
#' @param x A numeric vector of data
#' @param df1 The degrees of freedom of the numerator of the F distribution
#' @param df2 The degrees of freedom of the denominator of the F distribution
#'
#' @return The log-likelihood of \code{x} under the F distribution

#' @importFrom stats df
#' @export
log_likelihood_f = function(x, df1, df2) {
  return(sum(df(x, df1, df2, log = T)))
}


#' This function calculates the log-likelihood of a sample \code{x} under the t distribution, with degrees of freedom \code{df1} and \code{df2}.
#'
#' @param x A numeric vector of data
#' @param df The degrees of freedom of the t distribution
#'
#' @return The log-likelihood of \code{x} under the t distribution
#'
#' @examples
#' log_likelihood_t(rt(100, 5), 5)
#'
#' @importFrom stats dt
#' @export
log_likelihood_t = function(x, df) {
  return(sum(dt(x, df, log = T)))
}


#' This function calculates the sensitivity of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The sensitivity of the binary classification model
#'
#' @examples
#' sensitivity(c(1,1,0,1), c(1,0,0,1))
#'
#' @export
sensitivity = function(pred, truth) {
  true_positives <- sum((pred == 1) & (truth == 1))
  false_negatives <- sum((pred == 0) & (truth == 1))
  return(true_positives / (true_positives + false_negatives))
}


#' This function calculates the specificity of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The specificity of the binary classification model
#'
#' @examples
#' specificity(c(1,1,0,1), c(1,0,0,1))
#'
#' @export
specificity = function(pred, truth) {
  true_negatives <- sum((pred == 0) & (truth == 0))
  false_positives <- sum((pred == 1) & (truth == 0))
  return(true_negatives / (true_negatives + false_positives))
}

#' This function calculates the precision of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The precision of the binary classification model
#'
#' @examples
#' precision(c(1,1,0,1), c(1,0,0,1))
#'
#' @export
precision = function(pred, truth) {
  true_positives <- sum((pred == 1) & (truth == 1))
  false_positives <- sum((pred == 1) & (truth == 0))
  return(true_positives / (true_positives + false_positives))
}

#' This function calculates the recall of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The recall (sensitivity) of the binary classification model
#'
#' @examples
#' recall(c(1,1,0,1), c(1,0,0,1))
#'
#' @export
recall = function(pred, truth) {
  return(sensitivity(pred, truth))
}

#' This function calculates the accuracy of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The accuracy of the binary classification model
#'
#' @examples
#' accuracy(c(1,0,0,0), c(1,0,0,1))
#' @export
accuracy = function(pred, truth) {
  return(mean(pred == truth))
}


#' This function calculates the F1 score of a binary classification model based on the predicted labels \code{pred} and the true labels \code{truth}.
#'
#' @param pred A binary vector of predicted labels
#' @param truth A binary vector of true labels
#'
#' @return The F1 score of the binary classification model
#'
#' @examples
#' f1(c(1,1,0,1), c(1,0,0,1))
#'
#' @export
f1 = function(pred, truth) {
  prec <- precision(pred, truth)
  rec <- recall(pred, truth)
  return(2 * ((prec * rec) / (prec + rec)))
}


#' This function calculates the minimum sample size per group needed for a two-sample t-test, given the expected effect size \code{d} and the desired statistical power \code{power}.
#'
#' @param d The expected Cohen's d
#' @param power The desired statistical power (default 0.8)
#'
#' @return The minimum sample size per group needed for a two-sample t-test
#'
#' @examples
#' minimum_n_per_group(0.5)
#'
#' @importFrom pwr pwr.t.test
#' @export
minimum_n_per_group = function(d, power = 0.8) {
  require(pwr)
  return(pwr.t2n.test(d = d, power = power)$n)
}


#' This function calculates the R-squared statistic between a vector of predicted values \code{pred} and a vector of true values \code{truth}.
#'
#' @param pred A vector of predicted values
#' @param truth A vector of true values
#'
#' @return The R-squared statistic between \code{pred} and \code{truth}
#'
#' @examples
#' r2(c(1,2,3,4), c(1,3,2,5))
#'
#' @export
r2 = function(pred, truth) {
  return(1 - (sum((truth - pred)^2) / sum((truth - mean
                                           (truth))^2)))
}


#' This function calculates the adjusted R-squared statistic between a vector of predicted values \code{pred} and a vector of true values \code{truth}, given the number of model parameters \code{n_p}.
#'
#' @param pred A vector of predicted values
#' @param truth A vector of true values
#' @param n_p The number of model parameters, excluding the intercept
#'
#' @return The adjusted R-squared statistic between \code{pred} and \code{truth}
#'
#' @examples
#' adj_R2(c(1,2,3,4), c(1,3,2,5), 1)
#'
#' @export
adj_R2 = function(pred, truth, n_p) {
  n <- length(pred)
  r_squared <- r2(pred, truth)
  return(1 - ((1 - r_squared) * (n - 1) / (n - n_p - 1)))
}

