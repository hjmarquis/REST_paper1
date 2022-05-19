rm(list = objects())

# if(dir.exists("~/R-4.0.1/library"))
# {
#   .libPaths(c("~/R-4.0.1/library",.libPaths()))
#   np = 20
# }else{
#   np = 8
# }

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
source("source/function_benchmark.R")

load("data/analysis/paper1_elig_phecode.rda")
load("result/analysis/paper_surg_analysis_bench.rda")

# cl = makeCluster(getOption("cl.cores", np))
# registerDoParallel(cl)
# registerDoRNG(seed = 531)
# 
# if(dir.exists("~/R-4.0.1/library"))
# {
#   clusterEvalQ(cl, .libPaths(c("~/R-4.0.1/library",.libPaths())))
# }

n = length(surv)
nrct = length(surv.rct)
B = 1000


grp.all = factor(rep("2006-2024", length(trt)))
boot.fit.bench = vector("list", B)
for (b in 1:B) 
{
  print(b)
  start.time = Sys.time()
  Boot.rct = sample(1:nrct, nrct, replace = T)
  Boot.ehr = sample(1:n, n, replace = T)
  
  boot.fit.bench[[b]] = ate.ehr.phecode.bench(surv[Boot.ehr],trt[Boot.ehr],
                                                   Z.main[Boot.ehr,], Z.inter[Boot.ehr,], Z.inter.1[Boot.ehr,],  
                                              grp.all, year.grp,
                                    surv.rct[Boot.rct], trt.rct[Boot.rct], Z.rct[Boot.rct,],
                                    or.pen = fit.bench$or.pen, 
                                    or.lambda = fit.bench$or.lambda,
                                    ps.pen = fit.bench$ps.pen,  
                                    ps.lambda = fit.bench$ps.lambda,
                                    pos.exclude = fit.bench$pos.exclude, 
                                    dr.pen = fit.bench$dr.pen, 
                                    dr.lambda  = fit.bench$dr.lambda
  )
  run.time = Sys.time() - start.time
  
  print(run.time)
}


save( boot.fit.bench, 
     file = "result/analysis/paper_surg_analysis_bench_boot.rda")


