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
source("source/function_benchmark.R")

load("data/analysis/paper1_elig_phecode.rda")

grp.all = factor(rep("2006-2024", length(trt)))

fit.bench = ate.ehr.phecode.bench(surv,trt,Z.main, 
                            Z.inter,
                            Z.inter.1,
                            grp.all, year.grp, 
                            surv.rct, trt.rct, Z.rct
)


save(fit.bench, file = "result/analysis/paper_surg_analysis_bench.rda"))
