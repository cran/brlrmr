\name{il}
\alias{il}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
il
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This provides the estimates using IL method as described in the reference.
}
\usage{
il(formula, data, parameter = NULL, family = binomial, alpha = 0.05,
interaction = FALSE, k = NULL, na.action)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
as in \code{\link{lm}}. The missing values of response are NA.
}
\item{data}{
%%     ~~Describe \code{data} here~~
as in \code{\link{lm}}. The first column of data is binary missing response. The missing observations are denoted by NA. The rest of the columns are covariates or explanatory variables.
}
  \item{parameter}{
%%     ~~Describe \code{parameter} here~~
The starting values of the parameters as (\eqn{\beta, \alpha}) where \eqn{\beta} is the parameters of original model and \eqn{\alpha} is the for the missing data model.
}
  \item{family}{
  as in \code{\link{glm}}. \code{brlrmr} currently supports only the \code{"binomial"} family with links.
  }
  \item{alpha}{
%%     ~~Describe \code{tau} here~~
This is used for upper 100(1 - alpha)\% point of standard Normal distribution. The default is 1.96.
}
  \item{interaction}{
%%     ~~Describe \code{interaction} here~~
TRUE or FALSE, whether to consider interaction in the missing data model. Currenly only one intercation between response and covariates is supported. FALSE by default.
}
  \item{k}{
%%     ~~Describe \code{k} here~~
Which covariate has interaction with response. Takes integer values. User must assign a value if interaction = TRUE.
}
\item{na.action}{
as in \code{\link{lm}}. Always set to na.pass. Note that setting any other value to na.action will remove the NA's from response and hence will break the code as this package is only intended for missing response data.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{n }{number of observations.}
  \item{nmissing }{the number of missing observations.}
  \item{missing.proportion }{proportion of missing observations.}
  \item{beta.hat }{parameter estimate of logsitic regression of y on x using IL method.}
  \item{beta.se.hat }{standard error using IL method.}
  \item{z.value }{Wald Z value using IL method.}
  \item{p.value }{p value using IL method.}
  \item{significance.beta }{is indicator output whether regressors are significant using IL method, 1 if significant and 0 if not significant.}
  \item{LCL }{Lower Confidence Limits of 100(1 - alpha)\% Confidence Intervals.}
  \item{UCL }{Upper Confidence Limits of 100(1 - alpha)\% Confidence Intervals.}
  \item{alpha.hat }{parameter estimate due to missing model using IL.}
  \item{alpha.se.hat }{standard error of the them.}
  \item{z.value.alpha }{Wald Z value for them.}
  \item{p.value.alpha }{p values for them.}
}
\references{
%% ~put references to the literature/web site here ~
Ibrahim, J. G. and Lipsitz, S. R. (1996). Parameter estimation from incomplete data in binomial regression when the missing data mechanism is nonignorable. \emph{Biometrics}, 52:1071--1078.
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

%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

\examples{
\dontrun{
#############################################
########### Simulated Example ###############
#############################################
data(simulated.data)  # load simulated data

# parameter definition
beta0 <- 1
beta1 <- 1
beta2 <- 1
beta3 <- 1
beta4 <- 1

# parameter definition for missing indicator
alpha0 <- -1.1
alpha1 <- -1
alpha2 <- 1
alpha3 <- 1
alpha4 <- 1
alpha5 <- -1

parameter <- c(beta0, beta1, beta2, beta3, beta4,
               alpha0, alpha1, alpha2, alpha3, alpha4, alpha5)

il(y ~ x1 + x2 + x3 + x4, data = simulated.data, parameter,
family = binomial(link = "logit"), na.action = na.pass)
}

\dontrun{
#############################################
##### Real data example with separation #####
#############################################
data(nhanes)  # load nhanes data
il(hyp ~ age2 + age3, data = nhanes, family = binomial(link = "logit"), na.action = na.pass)
# IL method encounters seperation
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
