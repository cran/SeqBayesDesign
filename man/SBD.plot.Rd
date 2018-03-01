\name{SBD.plot}
\alias{SBD.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting the sequential design points
%%  ~~function to do ... ~~
}
\description{Plotting the sequential design points from the sequential Bayesian design.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SBD.plot(dat, Obj, cri = c("C-optimality", "D-optimality"),
         y.at = seq(0.3, 1, 0.1))
}

\arguments{
  \item{dat}{the information about the input data. See \code{\link{dat.setting}}.
}
  \item{Obj}{The output obtained from \code{\link{SBD.sim}}.
}
  \item{cri}{Selection of the criterion, "C-optimality" and "D-optimality".
}
  \item{y.at}{the plot range of y-axis.

  }
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
SBD.plot(dat.CFT, sim.res, cri = "C-optimality")
}
\donttest{
sim.res <- SBD.sim(12, dat.CFT$Nobs, dat.CFT, prior, "uniform", n.int=50000, n.sample=500, 
           init[[2]], "lnor", dat.CFT$mu.fun, candidate = seq(0.35, 0.75, 0.05), quantile=0.1)
SBD.plot(dat.CFT, sim.res, cri = "C-optimality")
}
}