\name{lifetime.model}
\alias{lifetime.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lifetime model determination
%%  ~~function to do ... ~~
}
\description{Model the data of a life testing. Provide the better model from the log-normal distribution or the Weibull distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lifetime.model(dat, initials)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{the information about the input data. See \code{\link{dat.setting}}.
  }
  \item{initials}{initial values of parameters for parameter estimation
  }
}

\value{
\item{lifetime.model}{the appropriate distribution.}
\item{est.par}{esitmates of model parameters.}
\item{log.likelihood}{the value of the log-likelihood.}
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.CFT <- dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL, 
                       Cen.time = 200000, mu.fun = "E-C", show = TRUE)
res <- lifetime.model(dat.CFT, c(0.01, 0.3, 0.7))
res
}

