\name{QICD.BIC}
\alias{QICD.BIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC for QICD on tuning parameter searching
}
\description{
Use the high dimensional BIC for quantile regression model on \code{QICD} algorithm, produces a plot and return a value for \code{lambda}
}
\usage{
QICD.BIC(y, x, beta = NULL, const = 6, tau, lambda, 
a = 3.7, funname = "scad", intercept = TRUE, 
thresh = 1e-06, maxin = 100, maxout = 20, 
plot.off = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response \code{y} as in \code{QICD}.
}
  \item{x}{
\code{x} matrix as in \code{QICD}.
}
  \item{beta}{
\code{beta} vector as in \code{QICD}
}
  \item{const}{
a parameter to adjust the high dimensional BIC. A positive numerical value.
}
  \item{tau}{
\code{tau} value as in \code{QICD}
}
  \item{lambda}{
a user supplied \code{lambda sequence}. A numerical vector, which will be used as a pool for tuning parameter searching 
}
  \item{a}{
\code{a} value as in \code{QICD}
}
  \item{funname}{
\code{funname} character vector as in \code{QICD}
}
  \item{intercept}{
\code{intercept} logical value as in \code{QICD}
}
  \item{thresh}{
\code{thresh} threshold as in \code{QICD}
}
  \item{maxin}{
\code{maxin} as in \code{QICD}
}
  \item{maxout}{
\code{maxout} as in \code{QICD}
}
  \item{plot.off}{
a logical value to control if a plot of QBIC vs. \code{lambda} will be produced. Default is FALSE and a plot will be given.
}
  \item{\dots}{
other argument that can be passed to \code{plot}
}
}
\details{
The function run \code{QICD} \code{nfolds} times. For each specific \code{lambda}, the QBIC will be produced for comparison. Claim that \code{cv.QICD} does NOT search for values of \code{a}. 
}
\value{
an object of class "BIC.QICD" is returned, which is a list with the components of QBIC.
%%  If it is a LIST, use
\item{lambda }{
the values of \code{lambda} used in the fits.
}
\item{HBIC}{
The high dimensional BIC is given-vector of length nlambda, as in \code{QICD} 
}
\item{nzero}{
number of non-zero coefficients at each \code{lambda}
}
\item{lambda.min}{
value of \code{lambda} that gives minimum \code{HBIC}.
}
%% ...
}
\references{
Peng,B and Wang,L. (2015)An Iterative Coordinate Descent Algorithm for High-dimensional
Nonconvex Penalized Quantile Regression, Journal of Computational and Graphical Statistics
\url{http://amstat.tandfonline.com/doi/abs/10.1080/10618600.2014.913516}
\doi{10.1080/10618600.2014.913516} 

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
\code{\link{QICD}},\code{\link{QICD.cv}}, \code{\link{QBIC}}
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
res.BIC=QICD.BIC(y, x, beta1,const=6, tau, 
lambda=seq(8,10,by=0.1), a,funname="scad",intercept=intercept)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QBIC }
\keyword{ High dimensional}
\keyword{Quantile regression}% __ONLY ONE__ keyword per line
