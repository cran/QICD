\name{allzero}
\alias{allzero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Judge an all zero vector
}
\description{
A logical function to judge if a vector is all zeros. If so, TRUE will be returned; otherwise, FALSE will 
be returned
}
\usage{
allzero(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a numerical vector
}
}
\details{
A simple tiny function embedded in
\code{QICD}\cr
}
\value{
TRUE will be returned if the vector contains all zeros, otherwise, FALSE will be returned
}

\author{
Bo Peng
}

\examples{
x=rep(0,10)
allzero(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ zeros}

