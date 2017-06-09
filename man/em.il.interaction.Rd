\name{em.il.interaction}
\alias{em.il.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
em.il.interaction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It is called by the main function fil and is for internal use.
}
\usage{
em.il.interaction(parameter, X, full.missing.data, observed.data, full.data, k, family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameter}{
%%     ~~Describe \code{parameter} here~~
The starting values of the parameters as (\eqn{\beta, \alpha}) where \eqn{\beta} is the parameters of original model and \eqn{\alpha} is the for the missing data model.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
The design matrix with the intercept column.
}
  \item{full.missing.data}{
%%     ~~Describe \code{full.missing.data} here~~
The augmented response with design matrix and missing indicator 1 for missing data.
}
  \item{observed.data}{
%%     ~~Describe \code{observed.data} here~~
The observed response with design matrix and missing indicator 0 for observed data.
}
  \item{full.data}{
%%     ~~Describe \code{full.data} here~~
The observed response, augmented response with corresponding design matrix and missing indicator 0 for observed data and 1 for missing data.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
If interaction is present in the missing data model, then the k is the column number of covariate matrix which has interaction with the response.
}
\item{family}{
as in \code{\link{glm}}. \code{brlrmr} currently supports only the  \code{"binomial"} family with links
}
}
\references{
%% ~put references to the literature/web site here ~
Bias Reduction in Logistic Regression with Missing Responses when the Missing Data Mechanism is Nonignorable.
}
