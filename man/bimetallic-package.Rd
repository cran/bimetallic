\name{bimetallic-package}
\alias{bimetallic-package}
\alias{bimetallic}
\docType{package}
\title{
Power for SNP analyses using silver standard cases
}
\description{
A power calculator for Genome-wide association studies (GWAs) with combined gold (error-free) 
and silver (erroneous) phenotyping per McDavid A, Crane PK, Newton KM, Crosslin DR, et al. (2011)
}
\details{
\tabular{ll}{
Package: \tab bimetallic\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-08-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
\code{\link{chisq.sim.factory}}, \code{\link{chisq.lytic.func}} and \code{\link{dlambda}}, \code{\link{is.legal}} all use common parameters.
chisq.sim.factory returns a function that repeatedly calculates chi-square statistics (X^2) and point estimates of allelic odds ratios (\code{\link{allelicOR}}) for a study specified by the simulation parameters.
chisq.lytic.func returns lambda, the non-centrality parameter of the asymptotic sampling distribution of (X^2), which happens to be non-central chi-square with two degrees of freedom.
dlambda returns the derivative of lambda with respect to \code{gam_ca} or \code{gam_co}.  If this is positive then power is increasing with the addition of a silver-standard subject.
is.legal determines if the parameters provided are \sQuote{legal} in the sense of inducing valid probabilities in affected and unaffected populations in the simulation
}
\author{
Andrew McDavid

Maintainer: Andrew McDavid <amcdavid@fhcrc.org>
}
\references{
See McDavid A, Crane PK, Newton KM, Crosslin DR, et al. (2011) for a full description of the model and parameters.
}
\keyword{ package }
\examples{
##Make a chisq simulator under a study scenario
sim = chisq.sim.factory(R = 4, gam_ca = 3, gam_co = 0,
	ppv = .8, npv = 1, homRR = 2.2, N_co = 1000,
	maf = .1, prev = .01, model = "mult")

##Run one realization of the simulation
sim()
}
