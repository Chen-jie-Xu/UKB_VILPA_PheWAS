# 03_Cox.R
source("00_data_preparation.R")
library(dplyr)
library(readr)
library(purrr)
library(survival)
library(broom)
library(forestploter)
library(grid)

# Outcomes selected from PheWAS
target_outcomes <- c(
  "Essential hypertension","Obesity","Chronic Kidney Disease, Stage III","GERD",
  "Major depressive disorder","Rheumatoid arthritis","Cancer of bronchus; lung",
  "Parkinson's disease","Chronic airway obstruction")

plot_order <- c(
  "Major depressive disorder","Parkinson's disease","Obesity",
  "Essential hypertension","Chronic Kidney Disease, Stage III","GERD",
  "Rheumatoid arthritis","Chronic airway obstruction","Cancer of bronchus; lung")

target_phe <- map %>%
  filter(phe_desc %in% target_outcomes) %>%
  distinct(phecode,phe_desc) %>%
  mutate(phe_desc=factor(phe_desc,levels=plot_order)) %>%
  arrange(phe_desc)
stopifnot(nrow(target_phe)==length(target_outcomes))
# VILPA tertiles
base_cox <- base %>%
  mutate( VILPA_tertile=cut(.data[[exposure_var]],c(-Inf,tertile_cuts,Inf),labels=c("T1","T2","T3")), VILPA_tertile=relevel(VILPA_tertile,ref="T1") )

fit_tertile <- function(phe,outcome){
  d <- make_surv_data(phe,base_cox)
  fit <- coxph(reformulate(c("VILPA_tertile",covariates), "Surv(followup,event)"),data=d)
  ne <- d %>%group_by(VILPA_tertile) %>% summarise(n=n(),events=sum(event),.groups="drop")
  hr <- tidy(fit,exponentiate=TRUE,conf.int=TRUE) %>%
    filter(grepl("^VILPA_tertile",term)) %>%
    transmute( VILPA_tertile=sub("VILPA_tertile","",term),HR=estimate,lower=conf.low,upper=conf.high,p=p.value)
  bind_rows(tibble(VILPA_tertile="T1",HR=1,lower=1,upper=1,p=NA_real_),hr
  ) %>%left_join(ne,by="VILPA_tertile") %>% mutate(phe_desc=outcome,phecode=phe,.before=1)
}

tertile_res <- map2_dfr(
  target_phe$phecode,as.character(target_phe$phe_desc),fit_tertile
) %>% mutate(
    p_fdr=ifelse(is.na(p),NA,p.adjust(p,method="fdr")),
    `HR (95% CI)`=ifelse(VILPA_tertile=="T1","1.00 (Reference)",sprintf("%.2f (%.2f, %.2f)",HR,lower,upper)),
    `Events/N`=paste0(events,"/",n) )

write_csv(tertile_res,file.path(out_dir,"Cox_tertile_results.csv"))

