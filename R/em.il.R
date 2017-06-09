#' @export
#'
em.il <-
function(parameter, X, full.missing.data, observed.data, full.data, family = family)
{  # parameter = (beta, alpha)

  p1 <- ncol(X)

  beta <- parameter[1:p1]
  alpha <- parameter[(p1 + 1):length(parameter)]


  linear.pred.y <- X %*% beta
  linear.pred.r <- cbind(X, full.missing.data[, 1]) %*% alpha


  likelihood.y <- exp(full.missing.data[, 1]*linear.pred.y)/(1 + exp(linear.pred.y))
  likelihood.r <- exp(full.missing.data[, (p1 + 1)]*linear.pred.r)/
    (1 + exp(linear.pred.r))
  likelihood.0 <- exp(0 * linear.pred.y)/(1 + exp(linear.pred.y))
  likelihood.1 <- exp(1 * linear.pred.y)/(1 + exp(linear.pred.y))


  weight.missing  <- (likelihood.y*likelihood.r)/
    (likelihood.0*likelihood.r + likelihood.1*likelihood.r)

  weight.observed <- rep(1, nrow(observed.data))


  weight <- c(weight.observed, weight.missing)


  full.x <- full.data[, 2:p1]
  full.y <- full.data[, 1]
  full.r <- full.data[, (p1 + 1)]



  glm.fit.y <- stats::glm(formula = full.y ~ full.x,
                   family = family,
                   weights = weight)
  # For Firth correction use brglm otherwise use glm2



  glm.fit.r <- stats::glm(formula = full.r ~ full.x + full.y,
                   family = family)
  # brglm is used to overcome the seperation problem


  beta.hat  <- stats::coef(glm.fit.y)
  alpha.hat <- stats::coef(glm.fit.r)
  Fisher    <- MASS::ginv(stats::vcov(glm.fit.y))  # Observed information matrix is inverse of
  # variance covariance matrix
  weights   <- glm.fit.y$weights
  Fisher.alpha <- MASS::ginv(stats::vcov(glm.fit.r))

  current.Q <- glm.fit.y$deviance/(-2) + glm.fit.r$deviance/(-2)

  parameter.hat <- c(beta.hat, alpha.hat)



  result <- list(Q = current.Q, parameter = parameter.hat,
                 Fisher = Fisher, weights = weights,
                 Fisher.alpha = Fisher.alpha)

  return(result)

}
