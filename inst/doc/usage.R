## ---- include = FALSE---------------------------------------------------------
library(historicalborrow)
library(dplyr)
library(posterior)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  fig.width = 7,
  fig.height = 5
)
set.seed(0)

## ---- paged.print = FALSE-----------------------------------------------------
library(historicalborrow)
library(dplyr)
set.seed(0)
data <- hb_sim_independent(
  n_continuous = 1,
  n_study = 3,
  n_group = 2,
  alpha = rep(1, 3),
  delta = 0.5,
  sigma = rep(1, 3),
  n_patient = 100
)$data %>%
  rename(
    outcome = response,
    trial = study,
    arm = group,
    subject = patient,
    factor1 = covariate_study1_continuous1,
    factor2 = covariate_study2_continuous1
  ) %>%
  mutate(
    trial = paste0("trial", trial),
    arm = paste0("arm", arm),
    subject = paste0("subject", subject)
  )
data

## -----------------------------------------------------------------------------
library(dplyr)
standardized_data <- hb_data(
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  covariates = c("factor1", "factor2")
)
standardized_data

## -----------------------------------------------------------------------------
distinct(
  standardized_data,
  study,
  study_label,
  group,
  group_label
) %>%
  select(
    study,
    study_label,
    group,
    group_label
  )

## -----------------------------------------------------------------------------
mcmc_pool <- hb_mcmc_pool(
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  # Can be continuous, categorical, or binary columns:
  covariates = c("factor1", "factor2"),
  # Raise these arguments for serious analyses:
  n_chains = 4,
  n_adapt = 2e3,
  n_warmup = 2e3,
  n_iterations = 4e3
)

mcmc_pool

## -----------------------------------------------------------------------------
mcmc_independent <- hb_mcmc_independent(
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  # Can be continuous, categorical, or binary columns:
  covariates = c("factor1", "factor2"),
  # Raise these arguments for serious analyses:
  n_chains = 4,
  n_adapt = 2e3,
  n_warmup = 2e3,
  n_iterations = 4e3
)

## -----------------------------------------------------------------------------
mcmc_hierarchical <- hb_mcmc_hierarchical(
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  # Can be continuous, categorical, or binary columns:
  covariates = c("factor1", "factor2"),
  # Raise these arguments for serious analyses:
  n_chains = 4,
  n_adapt = 2e3,
  n_warmup = 2e3,
  n_iterations = 4e3
)

## -----------------------------------------------------------------------------
hyperparameters <- hb_mcmc_mixture_hyperparameters(
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject"
)
hyperparameters

## -----------------------------------------------------------------------------
data_mixture <- dplyr::filter(data, trial == "trial3")
mcmc_mixture <- hb_mcmc_mixture(
  data = data_mixture, # only analyze current study
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  # Can be continuous, categorical, or binary columns:
  covariates = c("factor1", "factor2"),
  # Prior mixture components:
  m_omega = hyperparameters$m_omega,
  s_omega = hyperparameters$s_omega,
  p_omega = rep(1 / nrow(hyperparameters), nrow(hyperparameters)),
  # Raise these arguments for serious analyses:
  n_chains = 4,
  n_adapt = 2e3,
  n_warmup = 2e3,
  n_iterations = 4e3
)

## -----------------------------------------------------------------------------
hb_convergence(mcmc_hierarchical)

## -----------------------------------------------------------------------------
summary_hierarchical <- hb_summary(
  mcmc = mcmc_hierarchical,
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  covariates = c("factor1", "factor2"),
  eoi = c(0, 1),
  direction = c(">", "<")
)
summary_hierarchical

## -----------------------------------------------------------------------------
summary_pool <- hb_summary(
  mcmc = mcmc_pool,
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  covariates = c("factor1", "factor2")
)

summary_independent <- hb_summary(
  mcmc = mcmc_independent,
  data = data,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  covariates = c("factor1", "factor2")
)

hb_metrics(
  borrow = summary_hierarchical,
  pool = summary_pool,
  independent = summary_independent
)

## -----------------------------------------------------------------------------
summary_mixture <- hb_summary(
  mcmc = mcmc_mixture,
  data = data_mixture,
  response = "outcome",
  study = "trial",
  study_reference = "trial3",
  group = "arm",
  group_reference = "arm1",
  patient = "subject",
  covariates = c("factor1", "factor2")
)

hb_metrics(
  borrow = summary_mixture,
  pool = summary_pool,
  independent = summary_independent
)

## -----------------------------------------------------------------------------
hb_plot_borrow(
  borrow = summary_hierarchical,
  pool = summary_pool,
  independent = summary_independent
)

## -----------------------------------------------------------------------------
hb_plot_borrow(
  borrow = summary_mixture,
  pool = summary_pool,
  independent = summary_independent
)

## -----------------------------------------------------------------------------
hb_plot_group(
  borrow = summary_mixture,
  pool = summary_pool,
  independent = summary_independent
)

## -----------------------------------------------------------------------------
hb_plot_tau(mcmc_hierarchical)

