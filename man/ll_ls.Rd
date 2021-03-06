% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ls_functions.R
\name{ll_ls}
\alias{ll_ls}
\title{Log-likelihood}
\usage{
ll_ls(par, data, trans = FALSE)
}
\arguments{
\item{par}{vector c(piA, piB, muA1, muA0, muB1, sigma), c(piA, piB/(1-piA), muA1, muA0, muB1, sigma) if trans=TRUE.}

\item{data}{data frame containing columns y (positive outcome with zeros) and z (treatment).}

\item{trans}{boolean signifying if piB has been transformed.}
}
\value{
Log-likelihood for the latent stratification model.
}
\description{
Computes the log-likelihood for each of the four observational groups under the assumption of three strata and common variance.
}
\details{
For the input data frame, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
If z = 0, then the observation has not received treatment. \cr
Sometimes piB is transformed to relative proportions from absolute proportions. This transformation allows the reparameterization
of the piA and piB to allow constraint bounds between 0 and 1 in the optimization procedure. \cr
The log likelihoods are calculated based on equation 11-14 in the paper. Note that the equations presented in the paper are for
normal likelihood and for a single individual, thus corresponding adjustments have been made in the code to calculate the
log likelihood for the group. \cr
If you receive the warning \emph{Error in ll_ls(): Numeric overruns}, this means that the log likelihood has exceeded R's value storage
capacity, therefore storing it as infinity values. In this case, the value will be negative infinity, as likelihood does not exceed 1.
}
\examples{
sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
ll_ls(sim$par, sim$data)

# if the strata proportions are in relative sizes
ll_ls(sim$par, sim$data, trans=TRUE)

}
