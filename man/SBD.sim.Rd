\name{SBD.sim}
\alias{SBD.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation for the sequential Bayesian design
%%  ~~function to do ... ~~
}
\description{Given prior information and the true values of parameters,
provide the simulated new dataset by sequential Bayesian design.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SBD.sim(N.design, N.data, dat, prior, priorDis, n.int = 100000,
        n.sample = 1000, partrue, model, mu.fun=dat$mu.fun,
        candidate = seq(0.35, 0.75, 0.05), quantile=0.1, lambda = 0,
        L.lag = 100, transp = 0.5, show.acf = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N.design}{the number of next points should be determined.
  }
  \item{N.data}{number of the historical data.
  }
  \item{dat}{the information about the input data. See \code{\link{dat.setting}}.
  }
  \item{prior}{two parameters of a prior distribution for each model patameter.
  }
  \item{priorDis}{prior distribution. "uniform" is uniform distributed and "normal" is normal distributed.
  }
  \item{n.int}{the number of iteration in MCMC.
  }
  \item{n.sample}{the number of samples, which are used for evaluating the objective function by the Monte Carlo integration.
  }
  \item{partrue}{true values of parameter for generating a new observation in simulation.
  }
  \item{model}{the failure distribution. "lnor" is the log-normal distribution and "wei" is the Weibull distribution.
  }
  \item{mu.fun}{the relationship between location parameter and accelerating variable. See \code{\link{std.stress.fun}}.
  }
  \item{candidate}{The candidate points for the design.
  }
  \item{quantile}{The quantile of the lifetime.
  }
  \item{lambda}{Selection of the criterion. "0" is the C-optimality and "1" is the D-optimality.
  }
  \item{L.lag}{The largest lags of MCMC chain to be check the convergence. Default is 100.
  }
  \item{transp}{The model tranformation on the model parameters used in MCMC. Default is 0.5.
  }
  \item{show.acf}{logical; show the autocorrelation functions of the MCMC chains.
  }
}

\value{
\item{Bayes}{the Bayesian estimates of the historical data.}
\item{Final.opt}{The optimal design point and the values of objective functions.}
\item{data}{The augment data set by adding the next design point into the historical data.}
\item{hist.inf}{The determinant and the asymptotic variance of the historical data.}
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.CFT <- dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL, 
                       Cen.time = 200000, mu.fun = "E-C", show = TRUE)
init <- lifetime.model(dat.CFT, c(0.01, 0.3, 0.7))
prior <- c(0.000001, 0.1, 0.000001, 1, 4.5, 3)
# C-optimality
\dontshow{
sim.res <- SBD.sim(2, dat.CFT$Nobs, dat.CFT, prior, "uniform", n.int=10000, n.sample=50, 
           init[[2]], "lnor", dat.CFT$mu.fun, candidate = seq(0.35, 0.75, 0.5), quantile=0.1)
}
\donttest{
sim.res <- SBD.sim(12, dat.CFT$Nobs, dat.CFT, prior, "uniform", n.int=50000, n.sample=500, 
           init[[2]], "lnor", dat.CFT$mu.fun, candidate = seq(0.35, 0.75, 0.05), quantile=0.1)
}
}