

# RCT features: RCT and replicat
#----------------------------------------------
load("data/analysis/paper1_elig_phecode.rda")

colnames(Z.main)
colnames(Z.rct)

varname = c("Male", "Stage I", "Stage II", "Stage III",
            "Colon adhesion", "Right colon cancer","Left colon cancer", 
            "Sigmoid colon cancer")

Z.main = Z.main[,colnames(Z.rct)]
Z.main = cbind(Z.main,stage1=0,segRight=0)
Z.main[,"stage1"] = 1-Z.main[,"stage2"] 
Z.main[,"stage2"]  = Z.main[,"stage2"]  - Z.main[,"stage3"] 
Z.main[,"segRight"]  = 1-Z.main[,"segLeft"] -Z.main[,"segsigmoid"] 
Z.main = Z.main[,c("genderM", "stage1","stage2", "stage3", "adhesion", 
                   "segRight",  "segLeft",  "segsigmoid")]
or.ci = function(x,trt)
{
  fit = glm(x~trt, family = binomial)
  return(paste0(format(exp(coef(fit)[2]),digits=0,nsmall=2),
                " (", 
                paste(format(exp(confint(fit)[2,]),digits=0,nsmall=2),
                      collapse = ","),
                ')'))
}
n1 = sum(trt)
n0 =  sum(1-trt)
ehr1 = apply(Z.main*trt, 2, sum)
ehr0 = apply(Z.main*(1-trt), 2, sum)
ehr.tab = cbind(varname,
                paste(ehr1,' (',sapply(ehr1/n1*100,format,digits=0,nsmall=0),"%)",sep=''),
                paste(ehr0,' (',sapply(ehr0/n0*100,format,digits=0,nsmall=0),"%)",sep=''),
                apply(Z.main, 2, or.ci,trt=trt))

write.csv(ehr.tab, file = "report/paper surg/feature_rep.csv")


Z.rct = cbind(Z.rct,stage1=0,segRight=0)
Z.rct[,"stage1"] = 1-Z.rct[,"stage2"] 
Z.rct[,"stage2"]  = Z.rct[,"stage2"]  - Z.rct[,"stage3"] 
Z.rct[,"segRight"]  = 1-Z.rct[,"segLeft"] -Z.rct[,"segsigmoid"] 
Z.rct = Z.rct[,c("genderM", "stage1","stage2", "stage3", "adhesion", 
                  "segRight",  "segLeft",   "segsigmoid")]

n1 = sum(trt.rct)
n0 =  sum(1-trt.rct)
rct1 = apply(Z.rct*trt.rct, 2, sum)
rct0 = apply(Z.rct*(1-trt.rct), 2, sum)
rct.tab = cbind(varname,
                paste(rct1,' (',sapply(rct1/n1*100,format,digits=0,nsmall=0),"%)",sep=''),
                paste(rct0,' (',sapply(rct0/n0*100,format,digits=0,nsmall=0),"%)",sep=''),
                apply(Z.rct, 2, or.ci, trt=trt.rct))


write.csv(rct.tab, file = "report/paper surg/feature_rct.csv")
