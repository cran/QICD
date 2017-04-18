\name{checkloss}
\alias{checkloss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain the checkloss for quantile regression
}
\description{
A simple function. Return the checkloss for quantile regression.
}
\usage{
checkloss(res, tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
a numererical vector.
}
  \item{tau}{
quantile parameter in quantile regression
}
}
\value{
A checkloss value will be given, which is \deqn{\rho_{\tau}(u) = u(\tau-I(\tau < 0))}.
}
\references{
Peng,B and Wang,L. (2013)An Iterative Coordinate Descent Algorithm for High-dimensional
Nonconvex Penalized Quantile Regression, submitted
}
\author{
Bo Peng
}

\examples{
x=c(1,2,3)
checkloss(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ checkloss }
\keyword{ quantile regression}% __ONLY ONE__ keyword per line
