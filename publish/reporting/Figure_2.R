
load("data/analysis/paper1_elig_phecode.rda")
load("result/analysis/paper_surg_analysis_point.rda")
load("result/analysis/paper_surg_analysis_boot.rda")
local.dir = '' #  
rpdr.dir = paste0(local.dir,
                  "PHS data/Codified dictionary/rpdr_code_codebook_Marquis_i2b2.tsv")

rpdr.map = read.delim(rpdr.dir, stringsAsFactors = F)
rpdr.map = unique(rpdr.map[rpdr.map$feature_type == "DX",c("feature_id","bodysystem","feature_desc")])


var.sel = c("age","genderM", "raceW", "obesity", 
            "stage2", "stage3", "adhesion", 
            "segLeft",  "segsigmoid",
            "right.gen", "sigmoid.gen",
            "util_days", "util_yr",  
            "util_days_1yr", 
            "uniq_DX",
            "uniq_DX_1yr")
# Z.main = Z.main[,order(colnames(Z.main)%in% var.sel, decreasing = T)]
# Z.main = cbind(Z.main[,var.sel],
#                Z.main[, !(colnames(Z.main)%in% var.sel)])
varname = colnames(Z.main)
varname[match(var.sel, varname)] = c("Age", "Male", "White", "Obesity", 
                                     "Stage II", "Stage III", 
                                     "Colon adhesion", "Left colon cancer", 
                                     "Sigmoid colon cancer", 
                                     "General surgery code for right colon", 
                                     "General surgery code for sigmoid colon",
                                     "Healthcare Utilization, Overall", 
                                     "Healthcare Utilization Rate",
                                     "Healthcare Utilization, 1 Year",
                                     "Unique Diagnosis, Overall", "Unique Diagnosis, 1 Year")
varname = gsub("PheCode\\_", "PheCode:", varname)
pos.1yr = grep("1yr", varname)
pos.phecode = grep("PheCode", varname)
varname = gsub("\\_1yr", '', varname)
dx.descr = rpdr.map$feature_desc[match(varname[pos.phecode], rpdr.map$feature_id)]
for(na.dx in which(is.na(dx.descr)))
{
  dx.descr[na.dx] = paste0(unique(rpdr.map$bodysystem[grep(varname[pos.phecode[na.dx]],
                                                           rpdr.map$feature_id)]),
                           ": ", 
                           paste(unique(rpdr.map$feature_desc[grep(varname[pos.phecode[na.dx]],
                                                                   rpdr.map$feature_id)]),
                                 collapse = " or "))
}

varname[pos.phecode] = paste0(dx.descr, " (",
                              varname[pos.phecode], "), Overall")
varname[pos.1yr] = gsub("Overall", "1 Year", varname[pos.1yr])


coef.tab = data.frame(feature = rep(c("Time of Colectomy",varname), each = 2),
                       class = rep(c("Temporal","Demographic", "Cancer", "Healthcare", "Diagnosis"),
                                   c(2,8,14,10,length(pos.phecode)*2)),
                       arm = rep(c("OC", "LAC"), 1+ncol(Z.main)))

cox.tab = function(fit.death, Z.main, Z.inter)
{
  cox0609 = rep(fit.death$or.coef.off[1+1:ncol(Z.main)],
                each = 2)
  inter.main.pos = match(gsub("LAC_", '', colnames(Z.inter)), colnames(Z.main))
  cox0609[inter.main.pos*2] = cox0609[inter.main.pos*2] + fit.death$or.coef.off[1+ncol(Z.main)+1:ncol(Z.inter)]
  
  cox1417 =   cox1013 =cox0609
  
  base.var = c("year_0609","year_1013","year_1417","LAC_0609","LAC_1013","LAC_1417")
  for(j in which(! names(fit.death$or.coef) %in% base.var))
  {
    main.var = gsub("^LAC_|_0609$|_1013$|_1417$", '', names(fit.death$or.coef)[j])
    main.pos = which(colnames(Z.main) == main.var)
    if(grepl("^LAC_", names(fit.death$or.coef)[j]))
    {
      pos = 0
    }else{
      pos = -1:0
    }
    
    if(grepl("_0609$", names(fit.death$or.coef)[j]))
    {
      cox0609[main.pos*2 + pos] = cox0609[main.pos*2 + pos] + fit.death$or.coef[j]
    }else if(grepl("_1013$", names(fit.death$or.coef)[j]))
    {
      cox1013[main.pos*2 + pos] = cox1013[main.pos*2 + pos] + fit.death$or.coef[j]
    }else
    {
      cox1417[main.pos*2 + pos] = cox1417[main.pos*2 + pos] + fit.death$or.coef[j]
    }
  }
  pos.stage = which(colnames(Z.main)=="stage2")
  cox0609[pos.stage*2+1:2] = cox0609[pos.stage*2+1:2] + cox0609[pos.stage*2-1:0]
  cox1013[pos.stage*2+1:2] = cox1013[pos.stage*2+1:2] + cox1013[pos.stage*2-1:0]
  cox1417[pos.stage*2+1:2] = cox1417[pos.stage*2+1:2] + cox1417[pos.stage*2-1:0]
  cox0609 = c(fit.death$or.coef["year_0609"], 
              fit.death$or.coef["year_0609"] + fit.death$or.coef["LAC_0609"] + fit.death$or.coef.off[1], 
              cox0609)
  cox1013 = c(fit.death$or.coef["year_1013"], 
              fit.death$or.coef["year_1013"] + fit.death$or.coef["LAC_1013"] + fit.death$or.coef.off[1], 
              cox1013)
  cox1417 = c(fit.death$or.coef["year_1417"], 
              fit.death$or.coef["year_1417"] + fit.death$or.coef["LAC_1417"] + fit.death$or.coef.off[1], 
              cox1417)
  return(cbind(cox0609, cox1013, cox1417))
}

coef.tab = cbind(coef.tab, cox.tab(fit.death, Z.main, Z.inter))

boot.cox.tab = sapply(boot.fit.death, cox.tab, Z.main=Z.main, Z.inter=Z.inter)

cox.est = unlist(coef.tab[,c("cox0609", "cox1013", "cox1417")])
cox.sig = pmin(1,apply(sign(boot.cox.tab)!= sign(cox.est),1, mean) *2)
cox.sig[cox.est==0] = 1
coef.tab$cox0609.sig = cox.sig[1:nrow(coef.tab)]
coef.tab$cox1013.sig = cox.sig[1:nrow(coef.tab) + nrow(coef.tab)]
coef.tab$cox1417.sig = cox.sig[1:nrow(coef.tab) + 2*nrow(coef.tab)]

mean.overall = apply(array, margin, ...)

no.inter = which(apply((coef.tab[seq(1,nrow(coef.tab)-1, 2),c("cox0609", "cox1013", "cox1417")] - 
  coef.tab[seq(2,nrow(coef.tab), 2),c("cox0609", "cox1013", "cox1417")]) == 0,1,all))
coef.tab$arm[no.inter*2-1] = "Both"

Z.main[,"stage2"] = Z.main[,"stage2"] - Z.main[,"stage3"]
all.ave = c( length(year.grp)/3, apply(Z.main, 2, mean))
all.ave1 = c( sum(trt)/3, apply(Z.main[trt==1,], 2, mean))
all.ave0 = c( sum(1-trt)/3, apply(Z.main[trt==0,], 2, mean))

year.pos = year.grp=="2006-2009"
ave.0609 = c(sum(year.pos), apply(Z.main[year.pos,], 2, mean))/all.ave
ave1.0609 = c( sum(trt*year.pos)*2, apply(Z.main[trt==1 & year.pos,], 2, mean))/all.ave
ave0.0609 = c( sum((1-trt)*year.pos)*2, apply(Z.main[trt==0 & year.pos,], 2, mean))/all.ave

year.pos = year.grp=="2010-2013"
ave.1013 = c(sum(year.pos), apply(Z.main[year.pos,], 2, mean))/all.ave
ave1.1013 = c( sum(trt*year.pos)*2, apply(Z.main[trt==1 & year.pos,], 2, mean))/all.ave
ave0.1013 = c( sum((1-trt)*year.pos)*2, apply(Z.main[trt==0 & year.pos,], 2, mean))/all.ave

year.pos = year.grp=="2014-2017"
ave.1417 = c(sum(year.pos), apply(Z.main[year.pos,], 2, mean))/all.ave
ave1.1417 = c( sum(trt*year.pos)*2, apply(Z.main[trt==1 & year.pos,], 2, mean))/all.ave
ave0.1417 = c( sum((1-trt)*year.pos)*2, apply(Z.main[trt==0 & year.pos,], 2, mean))/all.ave

coef.tab$mean.0609 = c(rbind(ave0.0609, ave1.0609))
coef.tab$mean.0609[no.inter*2-1] = ave.0609[no.inter]
coef.tab$mean.1013 = c(rbind(ave0.1013, ave1.1013))
coef.tab$mean.1013[no.inter*2-1] = ave.1013[no.inter]
coef.tab$mean.1417 = c(rbind(ave0.1417, ave1.1417))
coef.tab$mean.1417[no.inter*2-1] = ave.1417[no.inter]

coef.tab = coef.tab[-no.inter*2,]
coef.tab = coef.tab[apply(coef.tab[,c("cox0609", "cox1013", "cox1417")]!=0| 
                            coef.tab$arm == "OC", 1, any),]

write.csv(coef.tab, file = "report/paper surg/visualization/cox_table_v3.csv",
          row.names = F)


coef.tab = data.frame(feature = c("Time of Colectomy",varname),
                      class = rep(c("Temporal","Demographic", "Cancer", "Healthcare", "Diagnosis"),
                                  c(1,4,7,5,length(pos.phecode))))


ps.tab = function(fit.death, Z.main, Z.time)
{
  ps1417 =   ps1013 =ps0609 = fit.death$ps.coef.off[-1]
  
  base.var = c("year_0609","year_1013","year_1417","LAC_0609","LAC_1013","LAC_1417")
  for(j in which(! colnames(Z.time) %in% base.var))
  {
    main.var = gsub("_0609$|_1013$|_1417$", '', colnames(Z.time)[j])
    main.pos = which(colnames(Z.main) == main.var)
    
    if(grepl("_0609$", colnames(Z.time)[j]))
    {
      ps0609[main.pos] = ps0609[main.pos] + fit.death$ps.coef[j+1]
    }else if(grepl("_1013$", colnames(Z.time)[j]))
    {
      ps1013[main.pos] = ps1013[main.pos] + fit.death$ps.coef[j+1]
    }else
    {
      ps1417[main.pos] = ps1417[main.pos] + fit.death$ps.coef[j+1]
    }
  }
  pos.stage = which(colnames(Z.main)=="stage2")
  ps0609[pos.stage+1] = ps0609[pos.stage+1] + ps0609[pos.stage]
  ps1013[pos.stage+1] = ps1013[pos.stage+1] + ps1013[pos.stage]
  ps1417[pos.stage+1] = ps1417[pos.stage+1] + ps1417[pos.stage]
  ps0609 = c(fit.death$ps.coef[2] + fit.death$ps.coef[1] + fit.death$ps.coef.off[1], 
              ps0609)
  ps1013 = c(fit.death$ps.coef[3] + fit.death$ps.coef[1] + fit.death$ps.coef.off[1], 
              ps1013)
  ps1417 = c(fit.death$ps.coef[4] + fit.death$ps.coef[1] + fit.death$ps.coef.off[1], 
              ps1417)
  return(cbind(ps0609, ps1013, ps1417))
}

coef.tab = cbind(coef.tab, ps.tab(fit.death, Z.main, Z.time))

boot.ps.tab = sapply(boot.fit.death, ps.tab, Z.main=Z.main, Z.time=Z.time)

ps.est = unlist(coef.tab[,c("ps0609", "ps1013", "ps1417")])
ps.sig = pmin(1,apply(sign(boot.ps.tab)!= sign(ps.est),1, mean) *2)
ps.sig[ps.est==0] = 1
coef.tab$ps0609.sig = ps.sig[1:nrow(coef.tab)]
coef.tab$ps1013.sig = ps.sig[1:nrow(coef.tab) + nrow(coef.tab)]
coef.tab$ps1417.sig = ps.sig[1:nrow(coef.tab) + 2*nrow(coef.tab)]

coef.tab$mean.0609 = ave.0609
coef.tab$mean.1013 = ave.1013
coef.tab$mean.1417 = ave.1417
coef.tab = coef.tab[apply(coef.tab[,c("ps0609", "ps1013", "ps1417")]!=0, 1, any) ,]


write.csv(coef.tab, file = "report/paper surg/visualization/ps_table.csv",
          row.names = F)
