\name{allelicOR}
\alias{allelicOR}
\title{
Calculate an odds ratio from a 2x3 contigency table of genotype X phenotype 
}
\description{
Find the odds ratio, per allele, in biallelic marker in a diploid individual.  The locus must be under Hardy-Weinburg at the margin, ie, alleles are independent, for this calculation to hold.
}
\usage{
allelicOR(xtable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xtable}{
A 2x3 contingency table of frequencies or counts of individuals.
}
}
\details{
The OR is calculated as (2*#AA.Case + #Aa.Case) * (2*#AA.Control + #Aa.Control) / 
((2*#aa.Case + #Aa.Case) * (2*#aa.Control + #Aa.Control))
}
\value{
The odds ratio, floating point.
}

\author{
Andrew McDavid
}

\seealso{
\code{\link{chisq.sim.factory}}, ~~~
}
\examples{
cont.table = matrix(c(5, 4, 1, 4, 3, 4), nrow=2, ncol=3, byrow=TRUE)
allelicOR(cont.table)
}

