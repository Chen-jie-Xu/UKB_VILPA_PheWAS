#1. PheWAS

## Cox-PheWAS:  Phecode v1.2

source("00_data_preparation.R")

library(survival)
library(broom)
library(EValue)


##  Cox-PheWAS

fml <- as.formula( paste( "Surv(followup, event) ~", paste(c(exposure_var, covariates), collapse = " + ")))

run_phe <- function(phe) {
  d <- make_surv_data(phe)
  n_event <- sum(d$event)
  if (n_event < min_events) return(NULL)
  fit <- tryCatch(
    coxph(fml, data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure_var) %>%
    transmute(
      phecode = phe,
      n_total = nrow(d),
      n_event = n_event,
      median_followup = median(d$followup),
      hr = estimate,
      conf.low,
      conf.high,
      p.value
    )
}

res <- map_dfr(sort(unique(phe_first$phecode)),run_phe)


## Add disease names + FDR 



## Phecode information
phe_info <- map %>%distinct(phecode, phe_desc, phe_group)

## E-value function
get_evalue <- function(hr, low, high) {
  x <- EValue::evalues.HR(est = hr, lo = low, hi = high, rare = TRUE)
  c(E_value_estimate = x["E-values", "point"],
    E_value_CI = if (hr >= 1) x["E-values", "lower"] else x["E-values", "upper"])}

## Calculate E-values
evalue <- t(mapply( get_evalue,res$hr,res$conf.low,res$conf.high))

## Final results
res <- res %>%
  left_join(phe_info, by = "phecode") %>%
  mutate(
    E_value_estimate = evalue[, "E_value_estimate"],
    E_value_CI = evalue[, "E_value_CI"],
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_by = p.adjust(p.value, method = "BY"),
    p_holm = p.adjust(p.value, method = "holm"),
    hr_ci = sprintf("%.2f (%.2f, %.2f)", hr, conf.low, conf.high),
    evalue_ci = sprintf("%.2f (%.2f)", E_value_estimate, E_value_CI)
  ) %>%
  arrange(p_fdr, p.value)

## 9. Save 
write_csv(res, "PheWAS_results.csv")
