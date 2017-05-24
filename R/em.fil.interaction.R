em.fil.interaction <-
function(parameter, X, full.missing.data, observed.data,
                               full.data, k, family = family)
{  # parameter = (beta, alpha)
  # X         = design matrix with intercepts

  p1 <- ncol(X)

  beta <- parameter[1:p1]
  alpha <- parameter[(p1 + 1):length(parameter)]


  linear.pred.y <- X %*% beta
  linear.pred.r <- cbind(X, X[, k + 1] * full.missing.data[, 1],
                         full.missing.data[, 1]) %*% alpha


  likelihood.y <- exp(full.missing.data[, 1]*linear.pred.y)/(1 + exp(linear.pred.y))
  likelihood.r <- exp(full.missing.data[, (p1 + 1)]*linear.pred.r)/(1 + exp(linear.pred.r))
  likelihood.0 <- exp(0 * linear.pred.y)/(1 + exp(linear.pred.y))
  likelihood.1 <- exp(1 * linear.pred.y)/(1 + exp(linear.pred.y))


  weight.missing  <- (likelihood.y*likelihood.r)/(likelihood.0*likelihood.r + likelihood.1*likelihood.r)

  weight.observed <- rep(1, nrow(observed.data))


  weight <- c(weight.observed, weight.missing)


  full.x1 <- full.data[, -c(1, p1 + 1)]
  full.y  <- full.data[, 1]
  full.r  <- full.data[, p1 + 1]


  brglm.fit.y <- brglm::brglm(formula = full.y ~ full.x1,
                       family = family,
                       weights = weight)
  # For Firth correction use brglm otherwise use glm2

  temp.y    <- full.y * full.x1[, k]
  glm.fit.r <- brglm::brglm(formula = full.r ~ full.x1  + temp.y + full.y,
                     family = family)
  # brglm is used to overcome the seperation problem


  alpha.hat       <- stats::coef(glm.fit.r)
  beta.hat.firth  <- stats::coef(brglm.fit.y)
  Fisher          <- brglm.fit.y$FisherInfo
  weights         <- brglm.fit.y$weights

  current.Q.firth <- brglm.fit.y$deviance/(-2) + glm.fit.r$deviance/(-2)


  parameter.hat.firth <- c(beta.hat.firth, alpha.hat)
  Fisher.alpha        <- glm.fit.r$FisherInfo

  result <- list(Q.firth = current.Q.firth, parameter.firth = parameter.hat.firth,
                 Fisher.firth = Fisher, weights = weights,
                 Fisher.firth.alpha = Fisher.alpha)

  return(result)

}
