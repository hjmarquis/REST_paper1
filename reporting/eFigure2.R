# ROC
dat.MAP = read.csv("result/phenotyping/MAP_CRC_ICDday.csv",
                   stringsAsFactors = F)

library(pROC)
library(ggplot2)

# All patients
roc.MAP = roc(dat.MAP$label,dat.MAP$MAP, direction = "<")
roc.ICD = roc(dat.MAP$label,dat.MAP$ICD, direction = "<")
roc.CUI = roc(dat.MAP$label,dat.MAP$CUI, direction = "<")

leg.lab = c(paste("MAP, AUC =", format(roc.MAP$auc,digits=0,nsmall=3)),
            paste("ICD, AUC =", format(roc.ICD$auc,digits=0,nsmall=3)),
            paste("NLP, AUC =", format(roc.CUI$auc,digits=0,nsmall=3)))
ggdf = data.frame( Specificity = c(rev(roc.MAP$specificities),
                                   rev(roc.ICD$specificities),
                                   rev(roc.CUI$specificities)),
                   Sensitivity = c(rev(roc.MAP$sensitivities),
                                   rev(roc.ICD$sensitivities),
                                   rev(roc.CUI$sensitivities)),
                   Method = rep(leg.lab,
                                c(length(roc.MAP$specificities),
                                  length(roc.ICD$specificities),
                                  length(roc.CUI$specificities)))
)

p.all = ggplot(ggdf, aes(x=Specificity, y = Sensitivity, color = Method))
p.all = p.all + geom_line() + scale_x_reverse()
p.all = p.all + labs(title = "All 171 Annotated Patients")
p.all = p.all + theme(legend.position="bottom",
                      legend.direction = "vertical")

# 
# png("figures/roc_MAP_CRC_ICDday_all.png")
# plot(roc.MAP, col = "red")
# lines(roc.ICD, col = "blue", lty = 2)
# lines(roc.CUI, col = "green", lty = 2)
# legend("bottomright", c(paste("MAP, AUC =", format(roc.MAP$auc,digits=0,nsmall=3)),
#                         paste("ICD, AUC =", format(roc.ICD$auc,digits=0,nsmall=3)),
#                         paste("NLP, AUC =", format(roc.CUI$auc,digits=0,nsmall=3))),
#        col = c("red","blue","green"),
#        lty = c(1,2,2))
# dev.off()

# Filter positive
dat.MAP.positive = dat.MAP[dat.MAP$ICD>0,]
roc.MAP = roc(dat.MAP.positive$label,dat.MAP.positive$MAP, direction = "<")
roc.ICD = roc(dat.MAP.positive$label,dat.MAP.positive$ICD, direction = "<")
roc.CUI = roc(dat.MAP.positive$label,dat.MAP.positive$CUI, direction = "<")

leg.lab = c(paste("MAP, AUC =", format(roc.MAP$auc,digits=0,nsmall=3)),
            paste("ICD, AUC =", format(roc.ICD$auc,digits=0,nsmall=3)),
            paste("NLP, AUC =", format(roc.CUI$auc,digits=0,nsmall=3)))
ggdf = data.frame( Specificity = c(rev(roc.MAP$specificities),
                                   rev(roc.ICD$specificities),
                                   rev(roc.CUI$specificities)),
                   Sensitivity = c(rev(roc.MAP$sensitivities),
                                   rev(roc.ICD$sensitivities),
                                   rev(roc.CUI$sensitivities)),
                   Method = rep(leg.lab,
                                c(length(roc.MAP$specificities),
                                  length(roc.ICD$specificities),
                                  length(roc.CUI$specificities)))
)

p.fp = ggplot(ggdf, aes(x=Specificity, y = Sensitivity, color = Method))
p.fp = p.fp + geom_line() + scale_x_reverse()
p.fp = p.fp + labs(title = "67 Annotated Patients with CRC ICD")
p.fp = p.fp + theme(legend.position="bottom",
                    legend.direction = "vertical")

# 
# png("figures/roc_MAP_CRC_ICDday_positive.png")
# plot(roc.MAP, col = "red")
# lines(roc.ICD, col = "blue", lty = 2)
# lines(roc.CUI, col = "green", lty = 2)
# legend("bottomright", c(paste("MAP, AUC =", format(roc.MAP$auc,digits=0,nsmall=3)),
#                         paste("ICD, AUC =", format(roc.ICD$auc,digits=0,nsmall=3)),
#                         paste("NLP, AUC =", format(roc.CUI$auc,digits=0,nsmall=3))),
#        col = c("red","blue","green"),
#        lty = c(1,2,2))
# dev.off()


library(gridExtra)
library(grid)

png("figures/paper surg/MAP_ROC.png", 
    width = 2200, height = 1400, res = 300)
grid.arrange(arrangeGrob(p.all, p.fp,
                         layout_matrix = matrix(1:2,1,2),
                         widths = c(1,1)
),
ncol = 1
)
dev.off()
