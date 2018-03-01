\name{dat.setting}
\alias{dat.setting}
\title{Information of historical data for sequential Bayesian designs
}
\description{Collect information of historical data for sequential Bayesian designs.
}
\usage{
dat.setting(dat, exp.type = c("ALT", "CFT.EC"), use.level, max.level = NULL,
            Cen.time, mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"), show = TRUE)
}
\arguments{
  \item{dat}{input of data.
}
  \item{exp.type}{experimental type. "ALT" is an accelerated life test and "CFT.EC" is a constant amplitude fatigue test.
}
  \item{use.level}{normal use condition of a accelerating variable.
}
  \item{max.level}{the maximum (highest) stress level. If \code{exp.type == "CFT.EC"}, \code{max.level} could be NULL.
}
  \item{Cen.time}{the censored time in life tests.
}
  \item{mu.fun}{the relationship between location parameter and accelerating variable. See \code{\link{std.stress.fun}}.
}
  \item{show}{logical; if TRUE, show the table of the data.
}
}
\value{
It collects all information of data, which is used for sequential Bayesian designs.
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL,
           Cen.time = 200000, mu.fun = "E-C", show = TRUE)
}
