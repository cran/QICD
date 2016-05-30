\name{QICD-package}
\alias{QICD-package}
\docType{package}
\title{
Estimation of coefficients of nonconvex penalized quantile regression model
}
\description{
Estimation of coefficients of nonconvex
penalized quantile regression model by using the Iterative Coordinate Descent (QICD)  algorithm.
This algorithm relies on a tuning parameter lambda that will be
chosen by both k-fold cross validation and high dimensional BIC for quantile regression model.
}
\details{
\tabular{ll}{
Package: \tab QICD\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-05-29\cr
License: \tab GPL-2\cr
}
This is a package to utilize the QICD algorithm on penalized quantile regression. Accepts \code{x,y, lambda} as predictor matrix, response variable and tuning parameter. Three main functions are included:
\code{QICD}\cr
\code{cv.QICD}\cr
\code{BIC.QICD}
for coefficients estimation and tuning parameter selection respectively.
Three other tiny functions are included as a supplement:
\code{allzero}\cr
\code{checkloss}\cr
\code{QBIC}\cr
}

\author{
Bo Peng

Maintainer: Bo Peng, peng0199@umn.edu
}
\references{
Peng,B and Wang,L. (2015)An Iterative Coordinate Descent Algorithm for High-dimensional
Nonconvex Penalized Quantile Regression. Journal of Computational and Graphical Statistics
\url{http://amstat.tandfonline.com/doi/abs/10.1080/10618600.2014.913516}
\doi{10.1080/10618600.2014.913516}

Lee, E. R., Noh, H. and Park. B. (2013) Model Selection via Bayesian Information
Criterion for Quantile Regression Models. Journal of the American Statistical Associa-
tion, preprint.
\url{http://www.tandfonline.com/doi/pdf/10.1080/01621459.2013.836975}
\doi{10.1080/01621459.2013.836975}


Wang,L., Kim, Y., and Li,R. (2013+)  Calibrating non-convex penalized regression in ultra-high dimension. To appear in Annals of Statistics.
\url{http://users.stat.umn.edu/~wangx346/research/nonconvex.pdf}

Fan, J. and Li, R.(2001) Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties. 
Journal of American Statistical Association, 1348-1360.
\url{http://orfe.princeton.edu/~jqfan/papers/01/penlike.pdf}

Zhang,C. (2010) Nearly Unbiase Variable Selection Under Minimax Concave Penalty. The Annals of Statistics, Vol. 38, No.2, 894-942
\url{http://arxiv.org/pdf/1002.4734.pdf}
}

\keyword{QICD}
\keyword{Majorization Minimization}
\keyword{Coordiante descent}
\keyword{QBIC}
\keyword{cross validation}
\examples{
x=matrix(rnorm(10000),50)
n=dim(x)[1]
p=dim(x)[2]
intercept=1
y=x[,1]+x[,7]+x[,9]+0.1*rnorm(n)
beta1=rep(0,p+intercept)
tau=0.5
a=2.7
res=QICD(y,x,beta1,tau,lambda=10,a,"scad",intercept=intercept)
}
