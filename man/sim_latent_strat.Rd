% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ls_functions.R
\name{sim_latent_strat}
\alias{sim_latent_strat}
\title{Simulate Data}
\usage{
sim_latent_strat(
  n = 1e+05,
  p = 0.5,
  piA = 0.05,
  piB = 0.1,
  muA1 = 5,
  muA0 = 4,
  muB1 = 5,
  sigma = 1
)
}
\arguments{
\item{n}{total sample size, the default value is 100000.}

\item{p}{treatment proportions, the default value is 0.5.}

\item{piA}{proportion of strata A, the default value is 0.05.}

\item{piB}{proportion of strata B, the default value is 0.10.}

\item{muA1}{mean for strata A that received treatment, Z = 1, the default value is 5.}

\item{muA0}{mean for strata A that did not received treatment, Z = 0, the default value is 4.}

\item{muB1}{mean for strata B that received treatment, Z = 1, the default value is 5.}

\item{sigma}{variance for all strata, the default value is 1.}
}
\value{
A list containing a data frame, a numeric value, and two vectors. \cr
The data frame \emph{data} contains the outcome variable y, treatment dummy z, and strata,
and the mean-centered effects-coded dummies for strata.\cr
The numeric value \emph{ATE} is the true average treatment effect. \cr
The vector \emph{par} consists of the true parameters of which the data is simulated. \cr
}
\description{
Simulates the data from the latent stratification model with four strata.
}
\details{
The four strata are defined as: \cr
A = positive under treatment and control, or always buyer. \cr
B = positive under treatment only, or influenced buyer. \cr
C = never positive, or never buyer. \cr
The model assumes that those who are positive under control, also known as defiers only doesn't exist.
The model also assumes that all strata share the same variance. \cr
The outcome y of strata A and B are generated through mixture models of normal distributions.
Strata A is generated using two normal distributions, one with mean muA1 and the other muA0.
Strata B is generated using a normal distribution with mean muB1 and 0, representing those in strata B won't be positive without treatment.
Strata C is 0 at all times. \cr
For the data frame in the output, column z is the dummy variable for treatment. If z = 1, then the observation has received treatment.
If z = 0, then the observation has not received treatment.
}
\examples{
sim = sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)
sim$par
# a vector piA 0.2 piB 0.1 muA1 5 muA0 4.5 mUB 3 sigma 0.3
sim$data
# a data frame containing outcome variable y, treatment dummy z, and strata. The first 5000 rows have z=1.
sim$ATE
# 4.0

}
