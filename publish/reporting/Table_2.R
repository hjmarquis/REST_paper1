load("result/analysis/paper_surg_analysis_point.rda")
load("result/analysis/paper_surg_analysis_boot.rda")
load("data/analysis/paper1_elig_phecode.rda")
est = "DR"
trt.level = c("LAC","Open")

pw.CI95 = function(x, name, est, time.grid)
{
  boot = sapply(x, function(x)
  {
    tmp = x
    for (i in 1:length(name))
    {
      tmp = tmp[[name[i]]]
    }
    stepfun(tmp$time, c(1, tmp[[est]]))(time.grid)
  })
  return(data.frame(time = time.grid, 
                    low = apply(boot, 1, quantile, prob = 0.025), 
                    up = apply(boot, 1, quantile, prob = 0.975)))
}

if(length(fit.death$pos.exclude)>0)
{
  yr.wgt = table(year.grp[-fit.death$pos.exclude])/length(year.grp[-fit.death$pos.exclude])
}else
{
  yr.wgt = table(year.grp)/length(year.grp)
}

fit.death$mu1.all$surv.ave = drop(sapply(fit.death$mu1.year, "[[", "DR") %*% yr.wgt)
fit.death$mu0.all$surv.ave = drop(sapply(fit.death$mu0.year, "[[", "DR") %*% yr.wgt) 

boot.mu1.all = sapply(boot.fit.death, function(x)  
  stepfun(x$mu1.all$time,c(1,drop(sapply(x$mu1.year, "[[", "DR") %*% yr.wgt)))(fit.death$mu1.all$time))
boot.mu0.all = sapply(boot.fit.death, function(x)  
  stepfun(x$mu0.all$time,c(1,drop(sapply(x$mu0.year, "[[", "DR") %*% yr.wgt)))(fit.death$mu0.all$time))


dim(boot.mu1.all)


tab.ehr = tab.diff = data.frame(time = c("All", levels(year.grp)),
                                OC = "", LAC = "", ATE = "",
                                stringsAsFactors = F)

format.cell = function(est,boot, digits = 3, type = c("emp","wald"))
{
  if(type[1]=="emp")
  {
    low = quantile(boot,prob=.025)
    up = quantile(boot,prob=.975)
  }else{
    se = sd(boot)
    low = est-1.96*se
    up = est+1.96*se
  }
  paste0(format(est, digits=0,nsmall = digits, trim = TRUE),
         " (", format(low, digits=0,nsmall = digits, trim = TRUE), 
         ',', format(up, digits=0,nsmall = digits, trim = TRUE), 
         ')')
}

full.oc = fit.death$mu0.all$surv.ave[which.min(abs(fit.death$mu0.all$time-5))]
full.lac = fit.death$mu1.all$surv.ave[which.min(abs(fit.death$mu1.all$time-5))]
full.ate = full.lac-full.oc

full.oc.boot = boot.mu0.all[which.min(abs(fit.death$mu0.all$time-5)),]
full.lac.boot = boot.mu1.all[which.min(abs(fit.death$mu1.all$time-5)),]
full.ate.boot = full.lac.boot - full.oc.boot

full.oc.diff = full.oc - fit.death$mu0.rct$surv[which.min(abs(fit.death$mu0.rct$time-5))]
full.lac.diff = full.lac - fit.death$mu1.rct$surv[which.min(abs(fit.death$mu1.rct$time-5))]
full.ate.diff = full.lac.diff - full.oc.diff

full.oc.diff.boot = full.oc.boot - sapply(boot.fit.death,function(x) 
  x$mu0.rct$surv[which.min(abs(x$mu0.rct$time-5))])
full.lac.diff.boot = full.lac.boot - sapply(boot.fit.death,function(x) 
  x$mu1.rct$surv[which.min(abs(x$mu1.rct$time-5))])
full.ate.diff.boot = full.oc.diff.boot - full.lac.diff.boot

tab.ehr$OC[1] = format.cell(full.oc, full.oc.boot)
tab.ehr$LAC[1] = format.cell(full.lac, full.lac.boot)
tab.ehr$ATE[1] = format.cell(full.ate, full.ate.boot)

tab.diff$OC[1] = format.cell(full.oc.diff, full.oc.diff.boot)
tab.diff$LAC[1] = format.cell(full.lac.diff, full.lac.diff.boot)
tab.diff$ATE[1] = format.cell(full.ate.diff, full.ate.diff.boot)

for (i in 1:3)
{
  yr.oc = fit.death$mu0.year[[i]]$DR[which.min(abs(fit.death$mu0.year[[i]]$time-5))]
  yr.lac = fit.death$mu1.year[[i]]$DR[which.min(abs(fit.death$mu1.year[[i]]$time-5))]
  yr.ate = yr.lac - yr.oc
  
  yr.oc.boot = sapply(boot.fit.death,function(x) 
    x$mu0.year[[i]]$DR[which.min(abs(x$mu0.year[[i]]$time-5))])
  yr.lac.boot = sapply(boot.fit.death,function(x) 
    x$mu1.year[[i]]$DR[which.min(abs(x$mu1.year[[i]]$time-5))])
  yr.ate.boot = yr.lac.boot - yr.oc.boot
  
  yr.oc.diff = yr.oc - fit.death$mu0.rct$surv[which.min(abs(fit.death$mu0.rct$time-5))]
  yr.lac.diff = yr.lac - fit.death$mu1.rct$surv[which.min(abs(fit.death$mu1.rct$time-5))]
  yr.ate.diff = yr.lac.diff - yr.oc.diff
  
  yr.oc.diff.boot = yr.oc.boot - sapply(boot.fit.death,function(x) 
    x$mu0.rct$surv[which.min(abs(x$mu0.rct$time-5))])
  yr.lac.diff.boot = yr.lac.boot - sapply(boot.fit.death,function(x) 
    x$mu1.rct$surv[which.min(abs(x$mu1.rct$time-5))])
  yr.ate.diff.boot = yr.oc.diff.boot - yr.lac.diff.boot
  
  tab.ehr$OC[i+1] = format.cell(yr.oc, yr.oc.boot)
  tab.ehr$LAC[i+1] = format.cell(yr.lac, yr.lac.boot)
  tab.ehr$ATE[i+1] = format.cell(yr.ate, yr.ate.boot)
  
  tab.diff$OC[i+1] = format.cell(yr.oc.diff, yr.oc.diff.boot)
  tab.diff$LAC[i+1] = format.cell(yr.lac.diff, yr.lac.diff.boot)
  tab.diff$ATE[i+1] = format.cell(yr.ate.diff, yr.ate.diff.boot)
}

write.csv(rbind(tab.ehr,tab.diff),
          file = "report/paper surg/result_table.csv",
          row.names = F)
