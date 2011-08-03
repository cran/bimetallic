\name{Power of Chi-square tests}
\alias{chisq.sim.factory}
\alias{is.legal}
\alias{chisq.lytic.func}
\alias{dlambda}
\title{
Power of Chi-square tests
}
\description{
Produce a function that returns a simulator of chi-square statistics and odds ratios, determine if
a set of parameters is legal, return the non-centrality parameter, lambda, of the chi-square distribution under a set of parameters, take a derivative of lambda with respect to number of silver-standard samples under a set of parameters.

}
\usage{
chisq.sim.factory(R, gam_ca, gam_co, ppv, npv, 
	homRR, N_co, maf, prev, model)
is.legal(R, gam_ca, gam_co, ppv, npv, 
	homRR, N_co, maf, prev, model)
chisq.lytic.func(R, gam_ca, gam_co, ppv, npv, 
	homRR, N_co, maf, prev, model)
dlambda(R, gam_ca, gam_co, ppv, npv, 
	homRR, N_co, maf, prev, model, diff="gam_ca")
}
\arguments{
  \item{R}{
Ratio of gold standard controls to gold standard cases
}
  \item{gam_ca}{
Ratio of silver standard cases to gold standard cases
}
  \item{gam_co}{
Ratio of silver standard controls to gold standard controls
}
  \item{ppv}{
The positive predictive value of the criteria used to identify silver standard cases, ie, P(Affected | Case)
}
  \item{npv}{
The negative predictive value of the criteria used to identify silver standard controls, ie, P(Unaffected | Control)
}
  \item{homRR}{
The relative risk of a homozygous genotype, ie, P(Affected | geno=AA) / P(Affected | geno=aa) 
}
  \item{N_co}{
The number of gold-standard controls
}
  \item{maf}{
The minor allele frequency of the risk locus
}
  \item{prev}{
Disease prevalence, ie, P(Affected)
}
  \item{model}{
Disease risk model, either one of \sQuote{dominant}, \sQuote{recessive}, \sQuote{multiplicative} or number giving the heterozygous relative risk.
}
  \item{diff}{Take the derivative of lambda with respect to \sQuote{gam_ca} (default) or \sQuote{gam_co}}
}
\details{
This function simulates the behavior of chi-square tests for independence in a genome-wide association study consisting of two cohorts.  The first cohort, the \dfn{gold standard} cohort, is assumed to have been classified into affected and unaffected without error.  The second cohort, the \dfn{silver standard} cohort, is assumed to have errors in disease classification subject to the \code{npv} and \code{ppv} arguments.  A distinction is drawn between \dfn{affected}, \dfn{unaffected} and \dfn{case} and \dfn{control}.  Affected and unaffected are the true disease status, which is observed in the gold standard cohort.  In the silver standard cohort the case/control criteria are observed instead, while the affected/unaffected status is latent.  

The numbers of gold and silver standard cases and controls are set by \code{N_co} and the ratios \code{R}, \code{gam_co} and \code{gam_ca}.  The genotype frequencies in the gold standard  cohort are governed by a genotypic disease  risk model and the parameters \code{homRR}, \code{maf}, \code{prev} and \code{model}.  The \code{model} argument may be a character vector of length one, either \sQuote{dominant}, \sQuote{recessive}, \sQuote{multiplicative} (or an unambigous abbreviation thereof), which links the heterozygous relative risk P(Affected | Aa) / P(Affected | aa)  to \code{homRR}, or it may be a floating point value directly specifying the heterozygous relative risk.

These functions offer a way to simulate power as well as calculate it asymptotically.   The distribution of chi-square statistics under the disease risk model, ppv, npv and case/control ratios is asymptotically distributed non-central chi-square.  \code{chisq.lytic.func} returns the non-centrality parameter associated with the model.  \code{dlambda} returns the first derivative of lambda.  Since power is increasing in lambda, if \code{dlambda} is positive, then an additional silver standard case (if \code{diff="gam_ca"}) increases power.
}
\value{
\code{chisq.sim.factory} returns a argument-less function that may be repeatedly called to sample from the simulated distribution.  This function returns a vector of length 2, consisting of the chi-square statistic and the point estimate of the allelic odds ratio.  Parameters are not checked for legality.

\code{is.legal} returns a boolean value indicating if the disease risk model parameters are legal, ie, induce genotype probabilities on [0,1], conditional on affected/unaffected status.

\code{chisq.lytic.func} returns \dfn{lambda}, the non-centrality parameter of the asymptotic sampling distribution of the chi-square test under the model.

\code{dlambda} returns the first derivative of lambda with respect to \code{gam_ca} or \code{gam_co}
}
\references{
McDavid A, Crane PK, Newton KM, Crosslin DR, et al. (2011) The Potential to Enhance Power of Genetic Association Studies through the Use of Silver Standard Cases Derived from Electronic Medical Records: Single Nucleotide Polymorphism Associations with Dementia in the Electronic Medical Records and Genomics (eMERGE) Network
}
\author{
Andrew McDavid
\email{amcdavidatfhcrcdotorg}
}

\seealso{
\code{\link{dlambda}}
}
\examples{
##Make a chisq simulator under a study scenario
sim = chisq.sim.factory(R = 4, gam_ca = 3, gam_co = 0,
	ppv = .8, npv = 1, homRR = 2.2, N_co = 1000,
	maf = .1, prev = .01, model = "mult")

##Run one realization of the simulation
sim()

##Run 100 realizations of the simulation
times = 100
chisq_or = as.data.frame(t(replicate(times, sim())))

##Find the number of times chisq_or$stat.X2[i] exceeded .01 significance
sig = .01
critval = qchisq(1-sig, 2)
nsucess = sum(chisq_or$stat > critval)

##Find the power
nsucess/times

##Compare to asymptotic
lambda = chisq.lytic.func(R = 4, gam_ca = 3, gam_co = 0,
	ppv = .8, npv = 1, homRR = 2.2, N_co = 1000,
	maf = .1, prev = .01, model = "mult")

1-pchisq(qchisq(1-sig, 2), 2, ncp=lambda)

## generate a multifactorial design
paramset = list(R = c(.5, 1, 2, 4), gam_ca = c(0, 1), gam_co = 0,
	ppv = c(.4, .6, .8), npv = .8, homRR = c(1.4, 3, 9), N_co = 1000,
	maf = .1, prev = .01, model = "mult")
paramgrid = do.call(expand.grid, c(paramset, stringsAsFactors=FALSE))

##call chisq.lytic and dlambda for each experiment
lambda = vector()
dl = vector()
for( i in 1:nrow(paramgrid)){
	lambda[i] = do.call(chisq.lytic.func, paramgrid[i,])
	dl[i] = do.call(dlambda, paramgrid[i,])
}

##bind it to the parameter data.frame 
##calculate the difference in lambda for gam_ca=0 vs gam_ca=1
paramgrid = cbind(paramgrid, lambda, dl)
param0 = subset(paramgrid, gam_ca==0)
param1 = subset(paramgrid, gam_ca==1)
all(paste(param0[,c(1, 3:10)])==paste(param1[,c(1, 3:10)]))
param0 = cbind(param0, finite_diff_dl = param1$lambda - param0$lambda)

##do they agree?
with(param0, cor(dl, finite_diff_dl))
}

