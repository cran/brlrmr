\name{em.fil.interaction}
\alias{em.fil.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
em.fil.interaction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It is called by the main function fil and is for internal use.
}
\usage{
em.fil.interaction(parameter, X, full.missing.data, observed.data, full.data, k, family)
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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
%% ~put references to the literature/web site here ~
Bias Reduction in Logistic Regression with Missing Responses when the Missing Data Mechanism is Nonignorable.
}
%\author{
%%  ~~who you are~~
%Arnab Kumar Maity
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
