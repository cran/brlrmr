#' @importFrom stats binomial qnorm contrasts is.empty.model model.matrix
#' model.response na.pass
#' @export
#'
fil <-
function(formula, data, parameter = NULL, family = binomial, alpha = 0.05, interaction = FALSE,
                k = NULL, na.action)
{

  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (family$family != "binomial")
    stop("families other than 'binomial' are not currently implemented")

  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "any")
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) {
    x <- NULL
  } else {
    x <- model.matrix(mt, mf, contrasts)
  }
  data <- cbind(y, x[, -1])

  tau <- qnorm(1 - alpha/2)
  n   <- nrow(data)  # number of observations
  p1  <- ncol(data)  # number of covariates + 1


  if(is.null(parameter))
    parameter <- rep(1, (2 * p1) + 1)

  original.data  <- data
  observed.data  <- cbind(subset(original.data, original.data[, 1] != "NA"), 0)
  missing.data   <- cbind(subset(original.data, is.na(original.data[, 1])), 1)
  nmissing       <- nrow(missing.data)  # no of missing observations



  # if y = 1 for the missing observations

  missing.data1 <- missing.data
  mi <- missing.data1[missing.data1[ , 1] <- 1]
  # replace all the reponses by 1




  # if y = 0 for the missing observations

  missing.data0 <- missing.data
  mi <- missing.data0[missing.data0[ , 1] <- 0]
  # replace all the reponses by 0




  full.data <- rbind(observed.data, missing.data1, missing.data0)
  full.missing.data <- rbind(missing.data1, missing.data0)

  x <- as.matrix(full.missing.data[, 2:p1])  # Covariates
  X <- cbind(1, x)       # Design Matrix
  p1 <- ncol(X)
  n.full <- nrow(x)


  if(interaction)
  {


    ############## Firth Correction ##################




    em.firth.result   <- em.fil.interaction(parameter, X, full.missing.data, observed.data,
                                            full.data, k, family = family)
    current.Q         <- em.firth.result$Q.firth
    current.parameter <- em.firth.result$parameter.firth


    error <- 1000
    count <- 1

    while(error > 0.0001)
    {

      em.firth.result   <- em.fil.interaction(parameter, X, full.missing.data, observed.data,
                                              full.data, k, family = family)
      next.Q            <- em.firth.result$Q.firth
      error             <- abs(next.Q - current.Q)
      current.Q         <- next.Q
      current.parameter <- em.firth.result$parameter.firth
      count             <- count + 1

    }


    beta.hat.firth <- current.parameter



    Fisher <- em.firth.result$Fisher.firth
    full.X <- cbind(1, full.data[, 2:p1])
    mu     <- boot::inv.logit(full.X %*% current.parameter[1:p1])
    w      <- em.firth.result$weights  # weights
    full.y <- full.data[, 1]
    full.r <- full.data[, (p1 + 1)]
    q      <- NULL
    W      <- diag(w)
    Ws     <- sqrt(W)
    H      <- Ws %*% full.X %*% chol2inv((chol(t(full.X) %*% W %*% full.X))) %*% t(full.X) %*% Ws


    for(j in 1:p1)
    {
      q[j] <- sum(w * (full.y - mu + H[j, j]*(0.5 - mu)) * full.X[, j])
    }

    s <- matrix(0, nrow = n.full, ncol = p1)

    for(i in 1:n.full)
    {
      for(j in 1:p1)
      {
        s[i, j] <- (full.y[i] - mu[i] + H[j, j]*(0.5 - mu[i])) * full.X[i, j]
      }  # S_i is p dimensional vector
    }

    second.term <- 0

    for(i in 1:n.full)
    {
      second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
    }

    Information <- Fisher - (second.term - q %*% t(q))
    se.matrix <- MASS::ginv(Information)


    beta.se.hat.firth <- sqrt(diag(se.matrix))


    Fisher.alpha <- em.firth.result$Fisher.firth.alpha
    full.X <- cbind(full.X, full.y * full.X[, k + 1], full.y)
    mu     <- boot::inv.logit(full.X %*% current.parameter[(p1 + 1):length(beta.hat.firth)])
    q      <- rep(0, p1 + 2)
    W  <- diag(w)
    Ws <- sqrt(W)
    H  <- Ws %*% full.X %*% chol2inv((chol(t(full.X) %*% W %*% full.X))) %*% t(full.X) %*% Ws


    for(j in 1:(p1 + 2))
    {
      q[j] <- sum(w * (full.r - mu + H[j, j]*(0.5 - mu)) * full.X[, j])
    }


    s <- matrix(0, nrow = n.full, ncol = (p1 + 2))

    for(i in 1:n.full)
    {
      for(j in 1:(p1 + 2))
      {
        s[i, j] <- (full.r[i] - mu[i] + H[j, j]*(0.5 - mu[i])) * full.X[i, j]
      }  # S_i is p dimensional vector
    }

    second.term <- 0

    for(i in 1:n.full)
    {
      second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
    }

    Information <- Fisher.alpha - (second.term - q %*% t(q))
    se.matrix <- MASS::ginv(Information)


    alpha.se.hat.firth <- sqrt(diag(se.matrix))

  } else {


    ############## Firth Correction ##################




    em.firth.result   <- em.fil(parameter, X, full.missing.data,
                                observed.data, full.data, family = family)
    current.Q         <- em.firth.result$Q.firth
    current.parameter <- em.firth.result$parameter.firth


    error <- 1000
    count <- 1

    while(error > 0.0001)
    {

      em.firth.result   <- em.fil(current.parameter, X, full.missing.data,
                                  observed.data, full.data, family = family)
      next.Q            <- em.firth.result$Q.firth
      error             <- abs(next.Q - current.Q)
      current.Q         <- next.Q
      current.parameter <- em.firth.result$parameter.firth
      count             <- count + 1

    }


    beta.hat.firth <- current.parameter



    Fisher <- em.firth.result$Fisher.firth
    full.X <- cbind(1, full.data[, 2:p1])
    mu     <- boot::inv.logit(full.X %*% current.parameter[1:p1])
    w      <- em.firth.result$weights  # weights
    full.y <- full.data[, 1]
    full.r <- full.data[, (p1 + 1)]
    q      <- NULL
    W      <- diag(w)
    Ws     <- sqrt(W)
    H      <- Ws %*% full.X %*% chol2inv((chol(t(full.X) %*% W %*% full.X))) %*% t(full.X) %*% Ws


    for(j in 1:p1)
    {
      q[j] <- sum(w * (full.y - mu + H[j, j]*(0.5 - mu)) * full.X[, j])
    }

    s <- matrix(0, nrow = n.full, ncol = p1)

    for(i in 1:n.full)
    {
      for(j in 1:p1)
      {
        s[i, j] <- (full.y[i] - mu[i] + H[j, j]*(0.5 - mu[i])) * full.X[i, j]
      }  # S_i is p dimensional vector
    }

    second.term <- 0

    for(i in 1:n.full)
    {
      second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
    }

    Information <- Fisher - (second.term - q %*% t(q))
    se.matrix <- MASS::ginv(Information)


    beta.se.hat.firth <- sqrt(diag(se.matrix))


    Fisher.alpha <- em.firth.result$Fisher.firth.alpha
    full.Xy <- cbind(full.X, full.y)
    mu <- boot::inv.logit(full.Xy %*% current.parameter[(p1 +1):length(beta.hat.firth)])
    q  <- rep(0, p1 + 1)
    W  <- diag(w)
    Ws <- sqrt(W)
    H  <- Ws %*% full.Xy %*% chol2inv((chol(t(full.Xy) %*% W %*% full.Xy))) %*% t(full.Xy) %*% Ws


    for(j in 1:(p1 + 1))
    {
      q[j] <- sum(w * (full.r - mu + H[j, j]*(0.5 - mu)) * full.Xy[, j])
    }

    s <- matrix(0, nrow = n.full, ncol = (p1 + 1))

    for(i in 1:n.full)
    {
      for(j in 1:(p1 + 1))
      {
        s[i, j] <- (full.r[i] - mu[i] + H[j, j]*(0.5 - mu[i])) * full.Xy[i, j]
      }  # S_i is p dimensional vector
    }

    second.term <- 0

    for(i in 1:n.full)
    {
      second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
    }

    Information <- Fisher.alpha - (second.term - q %*% t(q))
    se.matrix <- MASS::ginv(Information)


    alpha.se.hat.firth <- sqrt(diag(se.matrix))

  }


  return(list(n                  = n,
              nmissing           = nmissing,
              missing.proportion = nmissing/n,
              beta.hat           = beta.hat.firth[1:p1],
              beta.se.hat        = beta.se.hat.firth,
              z.value            = (abs(beta.hat.firth[1:p1] - 0)/beta.se.hat.firth),
              p.value            = 2 * (1 - stats::pnorm((abs(beta.hat.firth[1:p1] - 0)/beta.se.hat.firth))),
              significance.beta  = (abs(beta.hat.firth[1:p1] - 0)/beta.se.hat.firth) > tau,
              LCL                = beta.hat.firth[1:p1] - tau * beta.se.hat.firth,
              UCL                = beta.hat.firth[1:p1] + tau * beta.se.hat.firth,
              alpha.hat          = beta.hat.firth[(p1 +1):length(beta.hat.firth)],
              alpha.se.hat       = alpha.se.hat.firth,
              z.value.alpha      = (abs(beta.hat.firth[(p1 +1):length(beta.hat.firth)] - 0)/alpha.se.hat.firth),
              p.value.alpha      = 2 * (1 - stats::pnorm((abs(beta.hat.firth[(p1 +1):length(beta.hat.firth)] - 0)/alpha.se.hat.firth)))
  )
  )
}
