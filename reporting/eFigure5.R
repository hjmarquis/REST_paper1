rm(list = objects())

pak.list = c("glmnet", "ggplot2", "gridExtra", "grid", "survival", "cowplot")

for (pak in pak.list)
{
  yo = require(pak, character.only = T)
  if(!yo)
  {
    install.packages(pak,repos = "http://cran.us.r-project.org")
    require(pak, character.only = T)
  }
}
source("source/analysis/utility.R")
source("source/analysis/paper_surg_phecode.R")

load("result/analysis/paper_surg_analysis_point.rda")
load("result/analysis/paper_surg_analysis_boot.rda")
load("data/analysis/paper1_elig_phecode.rda")
load("result/analysis/paper_surg_analysis_bench.rda")
load("result/analysis/paper_surg_analysis_bench_boot.rda")

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

ci.mu1.rct = pw.CI95(boot.fit.death, "mu1.rct", "surv", fit.death$mu1.rct$time)
ci.mu0.rct = pw.CI95(boot.fit.death, "mu0.rct", "surv", fit.death$mu0.rct$time)
ci.mu1.bench = pw.CI95(boot.fit.bench, "mu1.all", "DR", fit.bench$mu1.all$time)
ci.mu0.bench = pw.CI95(boot.fit.bench, "mu0.all", "DR", fit.bench$mu0.all$time)

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

ci.mu1.ave = data.frame(time = fit.death$mu1.all$time, 
                        low = apply(boot.mu1.all, 1, quantile, prob = 0.025), 
                        up = apply(boot.mu1.all, 1, quantile, prob = 0.975))
ci.mu0.ave = data.frame(time = fit.death$mu0.all$time, 
                        low = apply(boot.mu0.all, 1, quantile, prob = 0.025), 
                        up = apply(boot.mu0.all, 1, quantile, prob = 0.975))

time.grid = sort(unique(surv[,1]))

fit.crude = list()

fit.crude$mu1.all = summary(survfit(surv~1, subset = trt==1),
                            time = time.grid)
fit.crude$mu0.all = summary(survfit(surv~1, subset = trt==0),
                            time = time.grid)

# Plot0: crude
#-------------------------------------------------------------

trt.level = c("LAC","Open")


p.title = paste0("Crude")
ggdf.all = data.frame(time = c(
                               0, fit.crude$mu1.all$time,
                               0, fit.crude$mu0.all$time),
                      surv = c(
                               1, fit.crude$mu1.all$surv,
                               1, fit.crude$mu0.all$surv),
                      low = c(
                              1, fit.crude$mu1.all$lower,
                              1, fit.crude$mu0.all$lower),
                      up = c(
                             1, fit.crude$mu1.all$upper,
                             1, fit.crude$mu0.all$upper),
                      trt = c(
                              rep(c(trt.level),c(length(fit.crude$mu1.all$time),
                                                 length(fit.crude$mu0.all$time))+1)),
                      stringsAsFactors = F)


tmp.p = ggplot(ggdf.all, aes(x = time, y = surv,color = trt))
tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt),
                            stat = "identity") 
tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", fill = "Surgical Type")
tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
# tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
#                                                      color = NA),
#                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
p.crude = tmp.p

# Plot1: benchmark
#-------------------------------------------------------------

p.title = paste0("Benchmark")
ggdf.bench = data.frame(time = c(
  0, fit.bench$mu1.all$time,
  0, fit.bench$mu0.all$time),
  surv = c(
    1, fit.bench$mu1.all$DR,
    1, fit.bench$mu0.all$DR),
  low = c(
    1, ci.mu1.bench$low,
    1, ci.mu0.bench$low),
  up = c(
    1, ci.mu1.bench$up,
    1, ci.mu0.bench$up),
  trt = c(
    rep(c(trt.level),c(length(fit.bench$mu1.all$time),
                       length(fit.bench$mu0.all$time))+1)),
  stringsAsFactors = F)


tmp.p = ggplot(ggdf.bench, aes(x = time, y = surv,color = trt))
tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt),
                            stat = "identity") 
tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", fill = "Surgical Type")
tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
# tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
#                                                      color = NA),
#                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
p.bench = tmp.p

# Plot2: benchmark
#-------------------------------------------------------------

p.title = paste0("Temporal-adjusted")
ggdf.death = data.frame(time = c(
  0, fit.death$mu1.all$time,
  0, fit.death$mu0.all$time),
  surv = c(
    1, fit.death$mu1.all$surv.ave,
    1, fit.death$mu0.all$surv.ave),
  low = c(
    1, ci.mu1.ave$low,
    1, ci.mu0.ave$low),
  up = c(
    1, ci.mu1.ave$up,
    1, ci.mu0.ave$up),
  trt = c(
    rep(c(trt.level),c(length(fit.death$mu1.all$time),
                       length(fit.death$mu0.all$time))+1)),
  stringsAsFactors = F)


tmp.p = ggplot(ggdf.death, aes(x = time, y = surv,color = trt))
tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt),
                            stat = "identity") 
tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", fill = "Surgical Type")
tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
# tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
#                                                      color = NA),
#                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
p.time = tmp.p

# Plot3: RCT
#-------------------------------------------------------------

p.title = paste0("RCT")
ggdf.rct = data.frame(time = c(0, fit.death$mu1.rct$time,
                               0, fit.death$mu0.rct$time),
                      surv = c(1, fit.death$mu1.rct$surv,
                               1, fit.death$mu0.rct$surv),
                      low = c(1, ci.mu1.rct$low,
                              1, ci.mu0.rct$low),
                      up = c(1, ci.mu1.rct$up,
                             1, ci.mu0.rct$up),
                      trt = rep(c(trt.level),c(length(fit.death$mu1.rct$time),
                                                 length(fit.death$mu0.rct$time))+1),
                      stringsAsFactors = F)
tmp.p = ggplot(ggdf.rct, aes(x = time, y = surv,color = trt))
tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt),
                            stat = "identity") 
tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", fill = "Surgical Type")
tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
# tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
#                                                      color = NA),
#                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
p.rct = tmp.p

g <- ggplotGrob(p.crude + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
png(paste0("figures/paper surg/analysis_CI_crude_bench_time_rct.png"), 
    width = 1200, height = 400)
grid.arrange(arrangeGrob(p.crude+ theme(legend.position="none"),
                         p.bench+ theme(legend.position="none"),
                         p.time+ theme(legend.position="none"),
                         p.rct+ theme(legend.position="none"),
                         layout_matrix = matrix(1:4,1,4),
                         widths = c(1,1,1,1)
),
legend,
ncol = 1,
heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()