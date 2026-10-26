# 05_KM.R
source("00_data_preparation.R")

library(survminer)
library(gridExtra)
library(grid)
library(scales)

# Outcomes selected from PheWAS
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
  distinct(phecode,phe_desc) %>%
  mutate(phe_desc=factor(phe_desc,levels=plot_order)) %>%
  arrange(phe_desc)

# Same tertile cut points as Cox
base_km <- base %>%mutate(tertile=cut(.data[[exposure_var]], c(-Inf,tertile_cuts,Inf),labels=c("T1","T2","T3")))

plot_km <- function(phe,outcome){
  d <- make_surv_data(phe,base_km)
  fit <- survfit(Surv(followup,event)~tertile,data=d)
  
  outcome <- recode(
    outcome,
    "Chronic Kidney Disease, Stage III"="Chronic kidney disease",
    "Cancer of bronchus; lung"="Lung cancer" )
  
  g <- ggsurvplot(
    fit,data=d,fun="event",
    risk.table=FALSE,conf.int=FALSE,censor=FALSE,
    xlim=c(0,8),break.time.by=2,
    palette=c("#F47c84","#7ec9f5","#f2de79"),
    legend.title="VILPA",legend.labs=c("T1","T2","T3"),
    xlab="Follow-up time (years)",
    ylab="Cumulative incidence (%)",
    title=outcome,
    ggtheme=theme_classic(base_size=12) )
  
  g$plot <- g$plot +
    scale_y_continuous(
      breaks=pretty_breaks(n=4),
      labels=label_number(scale=100,accuracy=.1)) +
    theme(
      plot.title=element_text(face="bold",size=12,hjust=.5),
      axis.title=element_text(face="bold",size=11),
      axis.text=element_text(size=10),
      legend.title=element_text(size=10),
      legend.text=element_text(size=9),
      legend.position="right"
    )
  
  g
}

km_list <- Map( plot_km,target_phe$phecode,as.character(target_phe$phe_desc))

combined_km <- arrangeGrob(grobs=lapply(km_list,`[[`,"plot"),ncol=3,nrow=3)

grid.newpage()
grid.draw(combined_km)

pdf(file.path(out_dir,"KM_curves.pdf"),width=12,height=10)
grid.draw(combined_km)
dev.off()