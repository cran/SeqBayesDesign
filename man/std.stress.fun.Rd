\name{std.stress.fun}
\alias{std.stress.fun}
\title{The value of location parameter in log-location and scale distributions
}
\description{The value of location parameter which is the funciton of the standardized accelerating variable with two parameters \eqn{A} and \eqn{B}.
}
\usage{
std.stress.fun(A, B, stlevels, mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"), dat)
}
\arguments{
  \item{A}{a value of parameter A.
  }
  \item{B}{a value of parameter B.
  }
  \item{stlevels}{
  vector of standardized stress levels.
  }
  \item{mu.fun}{the relationship between location parameter and accelerating variable. See "Details".
  }
  \item{dat}{the information about the input data. See \code{\link{dat.setting}}.
  }
}
\details{
Relationship "linear" is the linear function of stress levels (x). That is, \eqn{mu = A + B*x}.\cr
Relationship "Arrhenius" is the Arrhenius model. The standardized function is \eqn{mu = A + B/x}.\cr
Relationship "Inv.P" is the inverse-power law model. The standardized function is \eqn{mu = A + B*log(x)}.\cr
Relationship "E-C" is the S-N curve proposed by Epaarachchi and Clausen (2003).
}
\value{
The value of location parameter in log-location and scale distributions
}
\references{
Epaarachchi, J. A. and Clausen, P. D. (2003). "An empirical model for fatigue behavior
prediction of glass bre-reinforced plastic composites for various stress ratios and test
frequencies." \emph{Composites Part A: Applied science and manufacturing,} 34 (4), 313-326.
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.CFT <- dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL, 
                       Cen.time = 200000, mu.fun = "E-C", show = TRUE)
stlevel <- std.level(dat.CFT$use.level, dat.CFT$max.level, dat.CFT$x, mu.fun = "E-C")
res.lnorm <- par.MLE(dat.CFT, "lnor", "E-C", c(0.01, 0.3, 0.7))
std.stress.fun(res.lnorm$est.par[1], res.lnorm$est.par[2], stlevel, mu.fun = "E-C", dat.CFT)
}