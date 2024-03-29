% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_S4_generic_posterior.R
\name{posterior}
\alias{posterior}
\alias{posterior,RCModel-method}
\alias{posterior.RCModel}
\alias{posterior,Assessment-method}
\alias{posterior.Assessment}
\alias{RCMstan}
\title{Sample posterior of TMB models in SAMtool}
\usage{
posterior(x, ...)

\S4method{posterior}{RCModel}(
  x,
  priors_only = FALSE,
  laplace = FALSE,
  chains = 2,
  iter = 2000,
  warmup = floor(iter/2),
  thin = 5,
  seed = 34,
  init = "last.par.best",
  cores = chains,
  ...
)

\S4method{posterior}{Assessment}(x, priors_only = FALSE, ...)

RCMstan(RCModel, stanfit, sim, cores = 1, silent = FALSE)
}
\arguments{
\item{x}{An object of class \linkS4class{Assessment} or \linkS4class{RCModel}.}

\item{...}{Additional arguments to pass to \code{rstan::sampling} via \code{tmbstan::tmbstan}.}

\item{priors_only}{Logical, whether to set the likelihood to zero and sample the priors only.}

\item{laplace}{Logical, whether to do the Laplace approximation for random parameters.}

\item{chains}{The numer of MCMC chains.}

\item{iter}{The number of iterations for each chain, including warmup.}

\item{warmup}{The number of burnin iterations}

\item{thin}{The frequency at which iterations are kept (e.g., \code{5} saves every fifth iteration)}

\item{seed}{Seed for random number generator during the MCMC.}

\item{init}{The initial values of parameters for starting the MCMC chain. See \code{tmbstan::tmbstan}.}

\item{cores}{The number of cores for running in parallel, e.g., one core per MCMC chain. Used in \code{RCMstan}
for reconstructing the population.}

\item{RCModel}{An object of class \code{RCModel}}

\item{stanfit}{An object of class \code{stanfit} returned by \code{posterior}.}

\item{sim}{A matrix of \code{RCModel@OM@nsim} rows and 2 columns that specifies the samples used to update the
operating model. The first column specifies the chain and the second columns specifies the MCMC iteration.}

\item{silent}{Logical to indicate if progress messages should be printed to console.}
}
\value{
\code{posterior} returns an object of class \code{stanfit}. See \code{class?stanfit}.

\code{RCMstan} returns an updated \code{RCModel}.
}
\description{
A convenient wrapper function (\code{posterior}) to sample the posterior using MCMC in rstan
and returns a \code{stanfit} object for diagnostics. Use \code{RCMstan} to update the RCM and the enclosed operating model
with MCMC samples..
}
\section{Online Documentation}{

A vignette on the steps to run the MCMC is available on the openMSE
\href{https://openmse.com/tutorial-rcm/4-case-study-mcmc/}{website}.
}

\author{
Q. Huynh
}
