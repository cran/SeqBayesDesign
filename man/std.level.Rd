\name{std.level}
\alias{std.level}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardized stress level
}
\description{Transform the original scale of the accelerating variable to the standardized
stress levels.
}
\usage{
std.level(use.level, max.level, test.level, 
          mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{use.level}{normal use condition of a accelerating variable.
}
  \item{max.level}{the maximum (highest) stress level.
}
  \item{test.level}{the original scale of testing level.
}
  \item{mu.fun}{the relationship between location parameter and accelerating variable. See \code{\link{std.stress.fun}}.
}
}

\value{Return the standardized stress levels.
}

\seealso{
\code{\link{std.stress.fun}}.
}

\examples{
data(fiber)
data(use)
dat <- data.read(fiber, "CFT.EC")
dat.CFT <- dat.setting(dat, exp.type = "CFT.EC", use.level = use, max.level = NULL, 
                        Cen.time = 200000, mu.fun = "E-C", show = TRUE)
std.level(dat.CFT$use.level, dat.CFT$max.level, dat.CFT$x, mu.fun = "E-C")
}