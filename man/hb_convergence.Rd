% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hb_convergence.R
\name{hb_convergence}
\alias{hb_convergence}
\title{Check convergence diagnostics}
\usage{
hb_convergence(mcmc)
}
\arguments{
\item{mcmc}{A wide data frame of posterior samples returned by
\code{\link[=hb_mcmc_hierarchical]{hb_mcmc_hierarchical()}} or similar MCMC function.}
}
\value{
A data frame of summarized convergence diagnostics.
\code{max_rhat} is the maximum univariate Gelman/Rubin potential scale
reduction factor over all the parameters of the model,
\code{min_ess_bulk} is the minimum bulk effective sample size over the
parameters, and \code{min_ess_tail} is the minimum tail effective
sample size. \code{max_rhat} should be below 1.01, and the ESS metrics
should both be above 100 times the number of MCMC chains. If
any of these conditions are not true, the MCMC did not converge,
and it is recommended to try running the model for more saved
iterations (and if \code{max_rhat} is high, possibly more warmup
iterations).
}
\description{
Check the convergence diagnostics on a model.
}
\examples{
data <- hb_sim_pool(n_continuous = 2)$data
mcmc <- hb_mcmc_pool(
  data,
  n_chains = 1,
  n_adapt = 100,
  n_warmup = 200,
  n_iterations = 200
)
hb_convergence(mcmc)
}
\seealso{
Other mcmc: 
\code{\link{hb_mcmc_hierarchical}()},
\code{\link{hb_mcmc_independent}()},
\code{\link{hb_mcmc_mixture}()},
\code{\link{hb_mcmc_mixture_hyperparameters}()},
\code{\link{hb_mcmc_pool}()}
}
\concept{mcmc}
