rm(list = objects())

pak.list = c("glmnet", "ggplot2", "gridExtra", "grid", "survival"
             # ,"doParallel", "doRNG"
             )

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
load("result/analysis/paper_surg_analysis_point.rda")


n = length(surv)
nrct = length(surv.rct)
B = 1000

boot.fit = vector("list", B)
for (b in 1:B) 
{
  print(b)
  start.time = Sys.time()
  Boot.rct = sample(1:nrct, nrct, replace = T)
  Boot.ehr = sample(1:n, n, replace = T)
  
  boot.fit[[b]] = ate.ehr.phecode.by.year(surv[Boot.ehr],trt[Boot.ehr],
                                                   Z.main[Boot.ehr,], Z.inter[Boot.ehr,], Z.inter.1[Boot.ehr,], 
                                                   Z.time[Boot.ehr,], Z.inter.time[Boot.ehr,], Z.inter.1.time[Boot.ehr,], 
                                    year.grp[Boot.ehr],year.grp[Boot.ehr],
                                    surv.rct[Boot.rct], trt.rct[Boot.rct], Z.rct[Boot.rct,],
                                    or.pen.off = fit$or.pen.off, 
                                    or.lambda.off = fit$or.lambda.off, 
                                    or.lambda = fit$or.lambda,
                                    ps.pen.old.off = fit$ps.pen.old.off,   
                                    ps.lambda.old.off = fit$ps.lambda.old.off, 
                                    ps.lambda.old = fit$ps.lambda.old,
                                    ps.pen.off = fit$ps.pen.off,  
                                    ps.lambda.off = fit$ps.lambda.off, 
                                    ps.lambda = fit$ps.lambda,
                                    pos.exclude = fit$pos.exclude, 
                                    dr.pen = fit$dr.pen, 
                                    dr.lambda  = fit$dr.lambda
  )
  run.time = Sys.time() - start.time
  
  print(run.time)
}


save( boot.fit, 
     file = "result/analysis/paper_surg_analysis_boot.rda")


