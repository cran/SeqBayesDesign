\name{par.MLE}
\alias{par.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum likelihood estimates for lifetime model in accelerated life tests
}
\description{Maximum likelihood estimates for lifetime model in accelerated life tests (ALT).
}
\usage{
par.MLE(dat, model, mu.fun, initial, method = "BFGS", hessian = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{the information about the input data. See \code{\link{dat.setting}}.
%%     ~~Describe \code{dat} here~~
  }
  \item{model}{the failure distribution. "lnor" is the log-normal distribution and "wei" is the Weibull distribution.
%%     ~~Describe \code{model} here~~
  }
  \item{mu.fun}{the relationship between location parameter and accelerating variable. See \code{\link{std.stress.fun}}.
  }
  \item{initial}{initial values of parameter to be optimized over.
  }
  \item{method}{The method to be used. See \code{\link{optim}}.
  }
  \item{hessian}{Logical. Should a numerically differentiated Hessian matrix be returned?
  }
  \item{\dots}{additional arguments of \strong{optim()}. See \code{\link{optim}}.
  }
}

\value{
  \item{est.par}{the estimated values of parameters.
  }
  \item{log.likelihood}{the value of the log-likelihood function of the lifetime model
  }
}

\seealso{
  \code{\link{optim}}, \code{\link{lifetime.model}}
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.CFT <- dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL, 
                       Cen.time = 200000, mu.fun = "E-C", show = TRUE)
res.lnorm <- par.MLE(dat.CFT, "lnor", "E-C", c(0.01, 0.3, 0.7))
res.wei <- par.MLE(dat.CFT, "wei", "E-C", c(0.01, 0.3, 0.7))

res.lnorm
res.wei
}