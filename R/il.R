#' @importFrom stats binomial qnorm contrasts is.empty.model model.matrix
#' model.response na.pass
#' @export
#'
il <-
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

    sep <- 0  # separation indicator = 1 if separation, = 0 otherwise
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

      em.il.result      <- em.il.interaction(parameter, X, full.missing.data,
                                             observed.data, full.data, k,
                                             family = family)
      current.Q         <- em.il.result$Q
      current.parameter <- em.il.result$parameter


      error <- 1000
      count <- 1

      while(error > 0.0001)
      {

        em.il.result      <- em.il.interaction(current.parameter, X, full.missing.data,
                                               observed.data, full.data, k,
                                               family = family)
        next.Q            <- em.il.result$Q
        error             <- abs(next.Q - current.Q)
        current.Q         <- next.Q
        current.parameter <- em.il.result$parameter
        count             <- count + 1

      }


      beta.hat <- current.parameter

      if(sum(abs(beta.hat) > 10) > 0)
      {
        beta.hat = NA
        beta.se.hat = NA
        significance.beta = NA
        sep <- 1
      } else {
        w      <- em.il.result$weights  # weights
        n.full <- nrow(full.data)

        Fisher <- em.il.result$Fisher
        full.X <- cbind(1, full.data[, 2:p1])
        mu     <- boot::inv.logit(full.X %*% current.parameter[1:p1])
        full.y <- full.data[, 1]
        q      <- NULL

        for(j in 1:p1)
        {
          q[j] <- sum(w * (full.y - mu) * full.X[, j])
        }

        s <- matrix(0, nrow = n.full, ncol = p1)

        for(i in 1:n.full)
        {
          for(j in 1:p1)
          {
            s[i, j] <- (full.y[i] - mu[i]) * full.X[i, j]
          }  # S_i is p dimensional vector
        }

        second.term <- 0

        for(i in 1:n.full)
        {
          second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
        }

        Information <- Fisher - (second.term - q %*% t(q))
        se.matrix <- MASS::ginv(Information)


        beta.se.hat <- sqrt(diag(se.matrix))




        Fisher <- em.il.result$Fisher.alpha
        full.X <- cbind(full.X, full.y * full.X[, k + 1], full.y)
        mu     <- boot::inv.logit(full.X %*% current.parameter[(p1 + 1):length(beta.hat)])
        q      <- rep(0, p1 + 2)
        full.r <- full.data[, (p1 + 1)]



        for(j in 1:(p1 + 2))
        {
          q[j] <- sum(w * (full.y - mu) * full.X[, j])
        }



        s <- matrix(0, nrow = n.full, ncol = p1 + 2)

        for(i in 1:n.full)
        {
          for(j in 1:(p1 + 2))
          {
            s[i, j] <- (full.y[i] - mu[i]) * full.X[i, j]
          }  # S_i is p dimensional vector
        }

        second.term <- 0

        for(i in 1:n.full)
        {
          second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
        }

        Information <- Fisher - (second.term - q %*% t(q))
        se.matrix <- MASS::ginv(Information)


        alpha.se.hat <- sqrt(diag(se.matrix))
      }


    } else {

      em.il.result      <- em.il(parameter, X, full.missing.data,
                                 observed.data, full.data, family = family)
      current.Q         <- em.il.result$Q
      current.parameter <- em.il.result$parameter


      error <- 1000
      count <- 1

      while(error > 0.0001)
      {

        em.il.result      <- em.il(current.parameter, X, full.missing.data,
                                   observed.data, full.data, family = family)
        next.Q            <- em.il.result$Q
        error             <- abs(next.Q - current.Q)
        current.Q         <- next.Q
        current.parameter <- em.il.result$parameter
        count             <- count + 1

      }


      beta.hat <- current.parameter


      if(sum(abs(beta.hat) > 10) > 0)
      {
        beta.hat     <- NA
        beta.se.hat  <- NA
        alpha.se.hat <- NA
        alpha.se.hat <- NA
        sep          <- 1
      } else {

        w      <- em.il.result$weights  # weights
        n.full <- nrow(full.data)

        Fisher <- em.il.result$Fisher
        full.X <- cbind(1, full.data[, 2:p1])
        mu     <- boot::inv.logit(full.X %*% current.parameter[1:p1])
        full.y <- full.data[, 1]
        full.r <- full.data[, (p1 + 1)]
        q      <- NULL

        for(j in 1:p1)
        {
          q[j] <- sum(w * (full.y - mu) * full.X[, j])
        }

        s <- matrix(0, nrow = n.full, ncol = p1)

        for(i in 1:n.full)
        {
          for(j in 1:p1)
          {
            s[i, j] <- (full.y[i] - mu[i]) * full.X[i, j]
          }  # S_i is p dimensional vector
        }

        second.term <- 0

        for(i in 1:n.full)
        {
          second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
        }

        Information <- Fisher - (second.term - q %*% t(q))
        se.matrix <- MASS::ginv(Information)


        beta.se.hat <- sqrt(diag(se.matrix))

        Fisher.alpha <- em.il.result$Fisher.alpha
        q      <- rep(0, p1 + 1)


        for(j in 1:p1)
        {
          q[j] <- sum(w * (full.r - mu) * full.X[, j])
        }
        q[p1 + 1] <- sum(w * (full.r - mu) * full.y)

        s <- matrix(0, nrow = n.full, ncol = (p1 + 1))

        for(i in 1:n.full)
        {
          for(j in 1:p1)
          {
            s[i, j] <- (full.r[i] - mu[i]) * full.X[i, j]
          }  # S_i is p dimensional vector
          s[i, (p1 + 1)] <- (full.r[i] - mu[i]) * full.y[i]
        }

        second.term <- 0

        for(i in 1:n.full)
        {
          second.term <- second.term + w[i] * (s[i, ] %*% t(s[i, ]))
        }

        Information <- Fisher.alpha - (second.term - q %*% t(q))
        se.matrix <- MASS::ginv(Information)


        alpha.se.hat <- sqrt(diag(se.matrix))
      }


    }


    return(list(n                  = n,
                nmissing           = nmissing,
                missing.proportion = nmissing/n,
                beta.hat           = beta.hat[1:p1],
                beta.se.hat        = beta.se.hat,
                z.value            = (abs(beta.hat[1:p1] - 0)/beta.se.hat),
                p.value            = 2 * (1 - stats::pnorm((abs(beta.hat[1:p1] - 0)/beta.se.hat))),
                significance.beta  = (abs(beta.hat[1:p1] - 0)/beta.se.hat) > tau,
                LCL                = beta.hat[1:p1] - tau * beta.se.hat,
                UCL                = beta.hat[1:p1] + tau * beta.se.hat,
                alpha.hat          = beta.hat[(p1 +1):length(beta.hat)],
                alpha.se.hat       = alpha.se.hat,
                z.value.alpha      = (abs(beta.hat[(p1 +1):length(beta.hat)] - 0)/alpha.se.hat),
                p.value.alpha      = 2 * (1 - stats::pnorm((abs(beta.hat[(p1 +1):length(beta.hat)] - 0)/alpha.se.hat))),
                sep                = sep
    ))
  }
