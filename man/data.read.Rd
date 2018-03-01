\name{data.read}
\alias{data.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data Input
}
\description{Read a file in table format and creates a data frame from it.
}
\usage{
data.read(filename, exp.type = c("ALT", "CFT.EC"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{the name of the file which the data is to be read from. It could be the .csv file, \code{\link{read.csv}}. It could also be a data frame, with the necessary columns: \code{Stress}, \code{Lifetime}, \code{Censored}, and \code{Weight}. \code{Stress} is the testing stress levels, \code{Lifetime} is the lifetime of a unit, \code{Censored} is the indicator of a censored unit, and \code{Weight} is the observed weight.
}
  \item{exp.type}{experimental type. "ALT" is an accelerated life test and "CFT.EC" is a constant amplitude fatigue test. See "Details".
}
  \item{...}{additional arguments of \code{\link{read.csv}}.

  }
}
\details{
If \code{exp.type == "ALT"}, the order of the columns are \code{Stress}, \code{Lifetime},
 \code{Censored}, and \code{Weight}. \cr
If \code{exp.type == "CFT.EC"}, the order of the columns are \code{Stress}, \code{Lifetime}, \code{Censored}, \code{Weight}, \code{Rval}, \code{Frequency}, \code{Angle}, and \code{UTS}. 
\code{UTS} means the testing unit is failed at the first cycle if \code{UTS = 1}. All variables are needed for the S-N curve proposed by Epaarachchi and Clausen (2003).
}
\seealso{
  \code{\link{read.csv}}.
}

\references{
Epaarachchi, J. A. and Clausen, P. D. (2003). "An empirical model for fatigue behavior
prediction of glass bre-reinforced plastic composites for various stress ratios and test
frequencies." \emph{Composites Part A: Applied science and manufacturing,} 34 (4), 313-326.
}

\examples{
data(fiber)
data.read(fiber, exp.type = "CFT.EC")
}

