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



# Plot0: OR, IPW, DR OS up-to-5 years full cohort vs benchmark
#-------------------------------------------------------------

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

ci.mu1.rct = pw.CI95(boot.fit.death, "mu1.rct", "surv", fit.death$mu1.rct$time)
ci.mu0.rct = pw.CI95(boot.fit.death, "mu0.rct", "surv", fit.death$mu0.rct$time)
ci.mu1.all = pw.CI95(boot.fit.death, "mu1.all", est, fit.death$mu1.all$time)
ci.mu0.all = pw.CI95(boot.fit.death, "mu0.all", est, fit.death$mu0.all$time)
ci.mu1.ave = data.frame(time = fit.death$mu1.all$time, 
                        low = apply(boot.mu1.all, 1, quantile, prob = 0.025), 
                        up = apply(boot.mu1.all, 1, quantile, prob = 0.975))
ci.mu0.ave = data.frame(time = fit.death$mu0.all$time, 
                        low = apply(boot.mu0.all, 1, quantile, prob = 0.025), 
                        up = apply(boot.mu0.all, 1, quantile, prob = 0.975))
  
  p.title = paste0("All Years")
  ggdf.all = data.frame(time = c(0, fit.death$mu1.rct$time,
                             0, fit.death$mu0.rct$time,
                             0, fit.death$mu1.all$time,
                             0, fit.death$mu0.all$time),
                    surv = c(1, fit.death$mu1.rct$surv,
                             1, fit.death$mu0.rct$surv,
                             1, pmin(1,fit.death$mu1.all$surv.ave),
                             1, fit.death$mu0.all$surv.ave),
                    low = c(1, ci.mu1.rct$low,
                            1, ci.mu0.rct$low,
                            1, pmin(1,ci.mu1.ave$low),
                            1, ci.mu0.ave$low),
                    up = c(1, ci.mu1.rct$up,
                           1, ci.mu0.rct$up,
                           1, pmin(1,ci.mu1.ave$up),
                           1, ci.mu0.ave$up),
                    trt = c(rep(c(trt.level),c(length(fit.death$mu1.rct$time),
                                               length(fit.death$mu0.rct$time))+1),
                            rep(c(trt.level),c(length(fit.death$mu1.all$time),
                                               length(fit.death$mu0.all$time))+1)),
                    Source = c(rep(c("RCT", "EHR"),
                                   c(length(fit.death$mu1.rct$time)+
                                       length(fit.death$mu0.rct$time)+2,
                                     length(fit.death$mu1.all$time)+
                                       length(fit.death$mu0.all$time)+2))),
                    stringsAsFactors = F)
  
  ggdf.all$trt.alpha = 3-as.numeric(factor(ggdf.all$Source))
  
  tmp.p = ggplot(ggdf.all, aes(x = time, y = surv,  
                           linetype = Source, color = trt))
  tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt, 
                                     alpha = trt.alpha),
                                 stat = "identity") +  guides(alpha = FALSE) +
    scale_alpha(range = c(0.1, 0.4))
  tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
  tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", 
                        linetype = "Source", fill = "Surgical Type")
  tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
  # tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
  #                                                      color = NA),
  #                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
  p.all = tmp.p
  
  
  
  plist = ggdf.list = vector("list", 3)
  
  
  yr.levels = names(fit.death$mu1.year)
  
  for (i in 1:3) 
  { 
    
    ci.mu1 = pw.CI95(boot.fit.death, c("mu1.year",yr.levels[i]), est, 
                     fit.death$mu1.year[[yr.levels[i]]]$time)
    ci.mu0 = pw.CI95(boot.fit.death, c("mu0.year",yr.levels[i]), est, 
                     fit.death$mu0.year[[yr.levels[i]]]$time)
    
    p.title = paste0(yr.levels[i])
    ggdf = data.frame(time = c(0, fit.death$mu1.rct$time,
                               0, fit.death$mu0.rct$time,
                               0, fit.death$mu1.year[[i]]$time,
                               0, fit.death$mu0.year[[i]]$time),
                      surv = c(1, fit.death$mu1.rct$surv,
                               1, fit.death$mu0.rct$surv,
                               1, pmin(1,fit.death$mu1.year[[i]][[est]]),
                               1, fit.death$mu0.year[[i]][[est]]),
                      low = c(rep(NA,length(fit.death$mu1.rct$time)+
                                    length(fit.death$mu0.rct$time)+2),
                              1, pmin(1,ci.mu1$low),
                              1, ci.mu0$low),
                      up = c(rep(NA,length(fit.death$mu1.rct$time)+
                                   length(fit.death$mu0.rct$time)+2),
                             1, pmin(1,ci.mu1$up),
                             1, pmin(1,ci.mu0$up)),
                      trt = c(rep(c(trt.level),c(length(fit.death$mu1.rct$time),
                                                 length(fit.death$mu0.rct$time))+1),
                              rep(c(trt.level),c(length(fit.death$mu1.year[[i]]$time),
                                                 length(fit.death$mu0.year[[i]]$time))+1)),
                      Source = c(rep(c("RCT", "EHR"),
                                     c(length(fit.death$mu1.rct$time)+
                                         length(fit.death$mu0.rct$time)+2,
                                       length(fit.death$mu1.year[[i]]$time)+
                                         length(fit.death$mu0.year[[i]]$time)+2))),
                      stringsAsFactors = F)
    
    ggdf$trt.alpha = 3-as.numeric(factor(ggdf$Source))
    tmp.p = ggplot(ggdf, aes(x = time, y = surv,  
                                 linetype = Source, color = trt))
    tmp.p = tmp.p + geom_smooth(aes(ymin = low, ymax = up, fill = trt, 
                                    alpha = trt.alpha),
                                stat = "identity") +  guides(alpha = FALSE) +
      scale_alpha(range = c(0.1, 0.4))
    tmp.p = tmp.p + ylim(0.65, 1) + xlim(0,5)
    tmp.p = tmp.p +  labs(title = p.title,color = "Surgical Type", 
                          linetype = "Source", fill = "Surgical Type")
    tmp.p = tmp.p + ylab("Overall Survival") + xlab("Year(s) since Colectomy")
    # tmp.p = tmp.p + theme(plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
    #                                                      color = NA),
    #                       legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255)))
    
    plist[[i]] = tmp.p
    ggdf.list[[i]] = ggdf
    
    
  }

# The plot
g <- ggplotGrob(p.all + theme(legend.position="right",
                              legend.background = element_rect(fill=rgb(182/255, 182/255, 182/255))))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
lwidth <- sum(legend$widths)


pdf(paste0("figures/paper surg/analysis_CI_boot_poster.pdf"), 
    width = 15, height = 3,
    bg = rgb(182/255, 182/255, 182/255))
yo = cowplot::ggdraw(grid.arrange(arrangeGrob(p.all+ theme(legend.position="none",
                                                           plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
                                                                                          color = NA)),
                                              plist[[1]]+ theme(legend.position="none",
                                                                plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
                                                                                               color = NA)),
                                              plist[[2]]+ theme(legend.position="none",
                                                                plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
                                                                                               color = NA)),
                                              plist[[3]]+ theme(legend.position="none",
                                                                plot.background = element_rect(fill = rgb(182/255, 182/255, 182/255),
                                                                                               color = NA)),
                                              layout_matrix = matrix(1:4,1,4),
                                              widths = c(1,1,1,1)
),
legend,
nrow = 1,
widths = unit.c(unit(1, "npc") - lwidth, lwidth)
)) + theme(plot.background = element_rect(fill=rgb(182/255, 182/255, 182/255), color = NA))
dev.off()

g <- ggplotGrob(p.all + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
lheight <- sum(legend$height)
png(paste0("figures/paper surg/analysis_CI_boot.png"), 
    width = 1200, height = 400)
grid.arrange(arrangeGrob(p.all+ theme(legend.position="none"),
                         plist[[1]]+ theme(legend.position="none"),
                         plist[[2]]+ theme(legend.position="none"),
                         plist[[3]]+ theme(legend.position="none"),
                         layout_matrix = matrix(1:4,1,4),
                         widths = c(1,1,1,1)
),
legend,
ncol = 1,
heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()

pdf(paste0("figures/paper surg/analysis_CI_boot.pdf"), 
    width = 12, height = 4)
grid.arrange(arrangeGrob(p.all+ theme(legend.position="none"),
                         plist[[1]]+ theme(legend.position="none"),
                         plist[[2]]+ theme(legend.position="none"),
                         plist[[3]]+ theme(legend.position="none"),
                         layout_matrix = matrix(1:4,1,4),
                         widths = c(1,1,1,1)
),
legend,
ncol = 1,
heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()


setEPS()
postscript("figures/paper surg/analysis_CI_boot.eps",
           width = 12, height = 4)
grid.arrange(arrangeGrob(p.all+ theme(legend.position="none"),
                         plist[[1]]+ theme(legend.position="none"),
                         plist[[2]]+ theme(legend.position="none"),
                         plist[[3]]+ theme(legend.position="none"),
                         layout_matrix = matrix(1:4,1,4),
                         widths = c(1,1,1,1)
),
legend,
ncol = 1,
heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()

cairo_ps(file = "figures/paper surg/analysis_CI_boot.eps",
         width = 12, height = 4)
grid.arrange(arrangeGrob(p.all+ theme(legend.position="none"),
                         plist[[1]]+ theme(legend.position="none"),
                         plist[[2]]+ theme(legend.position="none"),
                         plist[[3]]+ theme(legend.position="none"),
                         layout_matrix = matrix(1:4,1,4),
                         widths = c(1,1,1,1)
),
legend,
ncol = 1,
heights = unit.c(unit(1, "npc") - lheight, lheight)
)
dev.off()


library(xlsx)

write.xlsx(ggdf.all, file = "report/paper surg/survival plots.xlsx", 
           sheetName = "All Years", row.names = FALSE)
write.xlsx(ggdf.list[[1]], file = "report/paper surg/survival plots.xlsx", 
           sheetName = "2006-2009", row.names = FALSE, append = TRUE)
write.xlsx(ggdf.list[[2]], file = "report/paper surg/survival plots.xlsx", 
           sheetName = "2010-2013", row.names = FALSE, append = TRUE)
write.xlsx(ggdf.list[[3]], file = "report/paper surg/survival plots.xlsx", 
           sheetName = "2014-2017", row.names = FALSE, append = TRUE)

