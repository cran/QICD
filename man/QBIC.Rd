\name{QBIC}
\alias{QBIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High dimensional BIC for quantile regression model
}
\description{
A high dimensional BIC will be returned specificall for quantile regression
}
\usage{
QBIC(y, X, beta, tau = 0.5, const = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response \code{y} as in \code{QICD}.
}
  \item{X}{
\code{x} matrix as in \code{QICD}.
}
  \item{beta}{
the coefficients vector for BIC calculation
}
  \item{tau}{
\code{tau} value as in \code{QICD}
}
  \item{const}{
a constant to adjust the BIC. A positive numerical value; default value is 6.
}
}
\details{
The high dimensional BIC for quantile regression model is 
\deqn{log(checkloss)+|S|log(log(n))C_n/n}
where S is the selected model in QICD, n is the number of obs, C_n is some positive constant which diverges to infinity as n increases. Actually, C_n is log(p)/\code{const}.
}

\value{
QBIC will be returned, which is a numerical value
}
\references{
Lee, E. R., Noh, H. and Park. B. (2013) Model Selection via Bayesian Information
Criterion for Quantile Regression Models. Journal of the American Statistical Associa-
tion, preprint.
\url{http://www.tandfonline.com/doi/pdf/10.1080/01621459.2013.836975}
\doi{10.1080/01621459.2013.836975}

Wang,L., Kim, Y., and Li,R. (2013+)  Calibrating non-convex penalized regression in ultra-high dimension. To appear in Annals of Statistics.
\url{http://users.stat.umn.edu/~wangx346/research/nonconvex.pdf}
}
\author{
Bo Peng
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{checkloss}}, \code{\link{QICD}}
}
\examples{
x=matrix(rnorm(1000),50)
n=dim(x)[1]
p=dim(x)[2]
intercept=1
y=x[,1]+x[,7]+x[,9]+0.1*rnorm(n)
beta1=rep(0,p+intercept)
tau=0.5
a=2.7
res=QICD(y,x,beta1,tau,lambda=10,a,"scad",intercept=intercept)
QBIC(y,cbind(x,rep(1,n)),res$beta_final,tau=tau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QBIC }
\keyword{ High dimensional }
\keyword{Quantile regression}% __ONLY ONE__ keyword per line
