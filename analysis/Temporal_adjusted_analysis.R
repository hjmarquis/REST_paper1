rm(list = objects())

if(dir.exists("~/R-4.0.1/library"))
{
  .libPaths(c("~/R-4.0.1/library",.libPaths()))
}

pak.list = c("glmnet", "ggplot2", "gridExtra", "grid", "survival")

for (pak in pak.list)
{
  yo = require(pak, character.only = T)
  if(!yo)
  {
    install.packages(pak,repos = "http://cran.us.r-project.org")
    require(pak, character.only = T)
  }
}
source("source/utility.R")
source("source/function_temporal_adjust.R")



load("data/analysis/paper1_elig_phecode.rda")

fit = ate.ehr.phecode.by.year(surv,trt,
                                       Z.main, Z.inter, Z.inter.1, 
                                       Z.time, Z.inter.time, Z.inter.1.time, 
                            year.grp,year.grp,
                            surv.rct, trt.rct, Z.rct

)


save( fit, file = "result/analysis/paper_surg_analysis_point.rda")

