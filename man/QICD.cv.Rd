\name{QICD.cv}
\alias{QICD.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for QICD 
}
\description{
Dose k-fold cross-validation for QICD, produces 
a plot and returns an appropriate tuning 
parameter \code{lambda}
}
\usage{
QICD.cv(y, x, beta = NULL, tau, lambda, nfolds = 10,
a = 3.7, funname = "scad", intercept = TRUE, thresh = 1e-06, 
maxin = 100, maxout = 20, mc.cores=getOption("mc.cores", 1L), 
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
  \item{tau}{
\code{tau} value as in \code{QICD}
}
  \item{lambda}{
a user supplied \code{lambda sequence}. A numerical vector, which will be used as a pool for tuning parameter searching 
}
  \item{nfolds}{
number of folds - default is 10.
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
  \item{mc.cores}{
  The number of cores to use for parallel computing, i.e. at most how many child processes will be run simultaneously. The option is initialized from environment variable MC_CORES if set. Must be at least one, and parallelization requires at least two cores.
}
  \item{plot.off}{
a logical value to control if a plot of prediction error vs. \code{lambda} will be produced. Default is FALSE and a plot will be given.
}
  \item{\dots}{
other argument that can be passed to \code{plot}
}
}
\details{
The function run \code{QICD} \code{nfolds} times. For each specific \code{lambda}, the average test prediction error will be produced for comparison. Claim that \code{QICD.cv} does NOT search for values of \code{a}.  
}
\value{
an object of class "cv.qicd" is returned, which is a list with the components of the cross-validation fit.
%%  If it is a LIST, use
\item{lambda }{
the values of \code{lambda} used in the fits.
}
\item{cvm}{
The mean cross-validated error-a vetor of length nlambda as in \code{QICD}
}
\item{cvsd}{
estimate of standard error of cvm.
}
\item{cvup}{
upper curve = cvm+cvsd.
}
\item{cvlo}{
upper curve = cvm-cvsd.
}
\item{nzero}{
number of non-zero coefficients at each \code{lambda}
}
\item{lambda.min}{
value of \code{lambda} that gives minimum \code{cvm}.
}
\item{lambda.1se}{
largest value of \code{lambda} such that error is within 1 standard error of the minimum.
}
%% ...
}
\references{
Peng,B and Wang,L. (2015)An Iterative Coordinate Descent Algorithm for High-dimensional
Nonconvex Penalized Quantile Regression. Journal of Computational and Graphical Statistics
\url{http://amstat.tandfonline.com/doi/abs/10.1080/10618600.2014.913516}
\doi{10.1080/10618600.2014.913516}
}
\author{
Bo Peng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QICD}},\code{\link{QICD.BIC}}
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
res.cv=QICD.cv(y, x, beta1, tau, lambda=seq(8,9,by=0.1),
nfolds=3, a, funname="scad",intercept=intercept)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cross Validation}
\keyword{Tuning parameter}% __ONLY ONE__ keyword per line
