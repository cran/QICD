\name{QICD}
\alias{QICD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a nonconvex penalized quantile regression model
}
\description{
Fit a nonconvex penalized quantile model via QICD algorithm. The estimation of the coefficients will be given. The regularization path is computed for the nonconvex penalties at a grid of values for the tuning parameter lambda. High dimentional BIC for quantile regression model (QBIC) and cross validation will be used as criterion on the tuning parameter searching.
}
\usage{
QICD(y, x, beta = NULL, tau, lambda, a = 3.7,
funname = "scad", intercept = TRUE, thresh = 1e-06, 
exclude = NULL, maxin = 100, maxout = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response variable. should be a numerical vector
}
  \item{x}{
input matrix, of dimension nobs x nvars; each row is an observation vector. In high dimensional cases, nvars could be larger than nobs
}
  \item{beta}{
initial value of coefficients. A numberical vector, with length of p if intercept is not included, and length of p+1 if intercept is included as the last entry, where p is the data dimention. The default value of beta is NULL, so an appropriate vector will be given automatically
}
  \item{tau}{
quantile parameter in quantile regression. A numerical value between 0 and 1.
}
  \item{lambda}{
a user supplied \code{lambda} sequence. A numerical vector. 
}
  \item{a}{
an unknown parameter in SCAD and MCD penalty functions. The default value is 3.7, suggested by Fan and Li (2001). 
}
  \item{funname}{
penalty function names. Character vector, which can be scad, mcp and lasso. 
}
  \item{intercept}{
Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE). Default value is TRUE.
}
  \item{thresh}{
Convergence threshold for QICD algorithm. Defaults value is 1E-6.
}
  \item{exclude}{
Indices of variables to be excluded from the model. Default is NULL. Can be logical vetor or numerical index specifying the exclued variables. 
}
  \item{maxin}{
Maximum number of inner coordiante descent iterations; default is 100.
}
  \item{maxout}{
Maximum number of outter Majoriaztion Minimization step (MM) iterations; default is 20.
}
}
\details{
The sequence of models indicated by \code{lambda} is fitted via the QICD algorithm. The QICD can substantially improves the computation speed in the p>>n setting. It combines the idea of the MM algorithm with that of the coordinate descent algorithm. For penalized quantile regression model, the objective function is 
\deqn{Q(\beta) = 1/n\sum_{i=1}^n\rho_{\tau}(Y_i-x_i^T\beta) +\sum_{j=1}^p p_{\lambda}(|\beta_j|)}
where \eqn{\rho_{\tau}} is the checkloss function for quantile regression. 
More specifically, we first replace the non-convex penalty function by its majorization function to create a surrogate objective function. Then we minimize the surrogate objective function with respect to a single parameter at each time and cycle through all parameters until convergence. For each univariate minimization problem, we only need
to compute a one-dimensional weighted median, which ensures fast computation.
}
\value{
An object with a list \code{QICD}.
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{beta_final }{
a matrix of coefficients. If intercept is included, diminsion is (p+1) x nlambda with the last row to be the intercepts, where nlamba is the length of \code{lambda} vector; if intercept is not included, diminsion is p x nlambda.
}
\item{lambda }{
the actual sequence of \code{lambda} values used
}
\item{df}{
The number of nonzero coefficients for each value of \code{lambda}
}
\item{dim}{
dimension of coefficient matrix (ices)
}
%% ...
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
\author{
Bo Peng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QICD.cv}}, \code{\link{QICD.BIC}}
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QICD}
\keyword{Coordiante Descent}
\keyword{Majorization Minimization}
\keyword{ nonconvex penalties }
\keyword{quantile regression}
\keyword{cross validation}
\keyword{HBIC}
% __ONLY ONE__ keyword per line
