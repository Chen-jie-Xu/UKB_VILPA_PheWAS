# 02_RCS.R

if (!exists("phe_first")) source("00_data_preparation.R")

library(rms)
library(ggplot2)
library(patchwork)

# Outcomes selected from PheWAS results
target_outcomes <- c(
  "Essential hypertension","Obesity","Chronic Kidney Disease, Stage III","GERD",
  "Major depressive disorder","Rheumatoid arthritis","Cancer of bronchus; lung",
  "Parkinson's disease","Chronic airway obstruction"
)

plot_order <- c(
  "Major depressive disorder","Parkinson's disease","Obesity",
  "Essential hypertension","Chronic Kidney Disease, Stage III","GERD",
  "Rheumatoid arthritis","Chronic airway obstruction","Cancer of bronchus; lung"
)

target_phe <- map %>%
  filter(phe_desc %in% target_outcomes) %>%
  distinct(phecode,phe_desc,phe_group) %>%
  mutate(phe_desc=factor(phe_desc,levels=plot_order)) %>%
  arrange(phe_desc)

fit_rcs <- function(phe,outcome){
  d <- make_surv_data(phe) %>% droplevels()
  if(sum(d$event)<min_events) return(NULL)
  
  knots <- quantile(d[[exposure_var]],c(.10,.50,.90),
                    na.rm=TRUE,names=FALSE)
  
  dd <- datadist(d %>% select(all_of(c(exposure_var,covariates))))
  dd$limits["Adjust to",exposure_var] <- 0
  assign("dd",dd,envir=.GlobalEnv)
  options(datadist="dd")
  
  spline <- paste0("rcs(",exposure_var,",c(",
                   paste(knots,collapse=","),"))")
  
  fml <- as.formula(paste(
    "Surv(followup,event) ~",
    paste(c(spline,covariates),collapse=" + ")
  ))
  
  fit <- cph(fml,data=d,x=TRUE,y=TRUE,surv=TRUE)
  
  pred <- do.call(
    Predict,c(list(object=fit),setNames(list(NA),exposure_var),
              list(ref.zero=TRUE,fun=exp))
  ) %>% as.data.frame()
  
  pred$phe_desc <- outcome
  pred$phecode <- phe
  
  a <- anova(fit)
  rn <- trimws(rownames(a))
  p_overall <- a[match(exposure_var,rn),"P"]
  p_nonlinear <- a[match("Nonlinear",rn),"P"]
  
  list(
    plot_data=pred,
    summary=tibble(
      phecode=phe,phe_desc=outcome,n=nrow(d),n_event=sum(d$event),
      p_overall=p_overall,p_nonlinear=p_nonlinear
    )
  )
}

rcs_list <- map2(target_phe$phecode,as.character(target_phe$phe_desc),fit_rcs) %>% compact()

plot_df <- map_dfr(rcs_list,"plot_data")
rcs_summary <- map_dfr(rcs_list,"summary")

# HR around median level(3 min/day)
hr3_df <- plot_df %>% group_by(phe_desc) %>%slice_min(abs(.data[[exposure_var]]-3),n=1,with_ties=FALSE) %>% ungroup()


# Plot
rename_outcome <- function(x) recode(
  x,
  "Chronic Kidney Disease, Stage III"="Chronic kidney disease",
  "Cancer of bronchus; lung"="Lung cancer"
)

plot_df$phe_desc <- rename_outcome(plot_df$phe_desc)
hr3_df$phe_desc <- rename_outcome(hr3_df$phe_desc)
rcs_summary$phe_desc <- rename_outcome(rcs_summary$phe_desc)

plot_order_rcs <- rename_outcome(plot_order)

plot_rcs <- function(outcome,show_x=FALSE,show_y=FALSE){
  pd <- filter(plot_df,phe_desc==outcome)
  hr <- filter(hr3_df,phe_desc==outcome)
  pv <- filter(rcs_summary,phe_desc==outcome)
  
  fp <- function(p) ifelse(p<.001,"<0.001",
                           paste0("= ",sprintf("%.2f",p)))
  
  ptxt <- paste0("P for overall ",fp(pv$p_overall),
                 "\nP for nonlinear ",fp(pv$p_nonlinear))
  
  hr$lab <- sprintf(
    "VILPA = %.1f\nHR = %.2f (%.2f, %.2f)",
    hr[[exposure_var]],hr$yhat,hr$lower,hr$upper
  )
  
  ggplot(pd,aes(.data[[exposure_var]],yhat)) +
    geom_ribbon(aes(ymin=lower,ymax=upper),  fill="#D8CBE8",alpha=.28) +
    geom_line(color="#6f4e9c",linewidth=1.05) +
    geom_hline(yintercept=1,linetype=2,color="grey40") +
    geom_point(data=hr,aes(.data[[exposure_var]],yhat), shape=21,size=2.7,fill="white",color="#C62828",inherit.aes=FALSE) +
    geom_label(data=hr,aes(.data[[exposure_var]],yhat,label=lab), hjust=-.08,vjust=1.35,size=3.4,fill="white", inherit.aes=FALSE) +
    annotate("text",x=Inf,y=Inf,label=ptxt, hjust=1.04,vjust=1.15,size=3.8) +
    labs(title=outcome, x=if(show_x) "VILPA (min/day)" else NULL, y=if(show_y) "Hazard ratio (95% CI)" else NULL) +
    scale_x_continuous(expand=expansion(mult=c(0,.02))) +
    scale_y_continuous(breaks=seq(0,1.2,.25), expand=expansion(mult=c(0,0))) +
    coord_cartesian(ylim=c(0,1.2),clip="off") +
    theme_classic(base_size=12)
}

plots <- map(seq_along(plot_order_rcs),\(i)
             plot_rcs(plot_order_rcs[i],
                      show_x=i>=7,
                      show_y=i %in% c(1,4,7)))

final_rcs <- wrap_plots(plots,ncol=3)
print(final_rcs)