\name{sparc}
\alias{sparc}
\title{
Training function of Semiparametric Generalized Linear Models
}
\description{
Pathwise Iterative Soft Thresholding Algorithm
}
\usage{
sparc(X, y, lambda = NULL, lambda.min.ratio=NULL, nlambda = NULL,
thol = 1e-4, max.ite = 1e4, alpha = sqrt(1/2))
}
\arguments{
\item{X}{
The \code{n} by \code{d} design matrix of the training set, where \code{n} is sample size and \code{d} is dimension.
}
\item{y}{
The \code{n}-dimensional response vector of the training set, where \code{n} is sample size.
}
\item{lambda}{
A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda. Supply instead a decreasing sequence of lambda values. sparc relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
}
\item{nlambda}{
The number of lambda values. The default value is 30.
}
\item{lambda.min.ratio}{
Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default is 0.1.
}
\item{thol}{
Stopping precision. The default value is 1e-4.
}
\item{max.ite}{
The number of maximum iterations. The default value is 1e4.
}
\item{alpha}{
The backtracking linea search parameter. The default value is sqrt(1/2).
}
}
\details{
We adopt pathwise Iterative Soft Thresholding Algorithm.
}
\value{
\item{w}{
The \code{d} by \code{nlambda} Regularization Path. The j-th column is the estimation coefficient corresponding to the j-th regularization parameter.
}
}
\references{
Y. Ning, Y. Chen, and H. Liu. "High Dimensional Semiparametric Generalized Linear Models", \emph{Technical Report}, 2013.
}
\author{
Tuo Zhao and Han Liu\cr
Maintainers: Tuo Zhao<tourzhao@gmail.com>
}
\seealso{
\code{\link{sparc-package}}
}
\examples{
## generating training data
n = 100
d = 200
set.seed(3)
X = matrix(rnorm(n*d),n,d)
y = 3*X[,1]+2*X[,2] + 1.5*X[,4] + rnorm(n)

## estimating models
out = sparc(X,y)
}