library(readxl)

dat.nlp = read.csv("data/EHR/Merck_PatientLeveldata_07062020.csv",
                     header = T)
# dat.note = read.delim("data/EHR/ColorectalCa_res_count.txt", sep = "|",
#                       header = F)
dat.all.ICD = read.csv("data/EHR/Merck_DIagnosis.csv",
                       stringsAsFactors = F, header = F,
                       fileEncoding = "UTF-8-BOM")
dat.label = read.csv("data/label/label200.csv",stringsAsFactors = F,
                     col.names = c("patientNum","Label","Comment"))[,-3]
dat.ICD = read_xlsx("data/EHR/Merck_AllPatients_Colorectal_Cancer_CodeCountv3.xlsx", sheet = "patients_ColonCancer_codeCount")
new.filter =   unlist(read_xlsx("dictionary/Merck_REST_new_colorectal_cancer_filter.xlsx", sheet = "New_filter_malignant"))

# Get all days with ICD
dat.all.ICD$Date = sapply(strsplit(dat.all.ICD$V2," "), function(x) x[1])
names(dat.all.ICD)[1]="patientNum"
dat.all.ICD = dat.all.ICD[,c("patientNum","Date")]
dat.all.ICD = unique(dat.all.ICD)


dat.MAP = data.frame(patientNum = unique(dat.all.ICD$patientNum),
                     note = c(table(dat.all.ICD$patientNum)),
                     ICD = 0,
                     CUI = 0)

dat.nlp = dat.nlp[is.element(dat.nlp$patientNum, dat.MAP$patientNum),]
dat.MAP$CUI[match(dat.nlp$patientNum, dat.MAP$patientNum)] = (dat.nlp$C0699790Y_CARCINOMA_OF_COLON+
  dat.nlp$C0007102Y_Malignant_tumor_of_colon +
  dat.nlp$C0007113Y_CARCINOMA_OF_RECTUM +
  dat.nlp$C0009402Y_CRC)

dat.label = dat.label[is.element(dat.label$patientNum,dat.MAP$patientNum),]
dat.MAP$label = NA
dat.MAP$label[match(dat.label$patientNum,dat.MAP$patientNum)] = dat.label$Label!=0

dat.ICD = dat.ICD[is.element(dat.ICD$code,new.filter),]
dat.ICD = dat.ICD[is.element(dat.ICD$patientNum,dat.MAP$patientNum),]
ICD.patientNum = unique(dat.ICD$patientNum)
ICD.count = diff(c(0,cumsum(dat.ICD$Total_colorectal_ICDCount)[cumsum(table(dat.ICD$patientNum))]))
dat.MAP$ICD[match(ICD.patientNum,dat.MAP$patientNum)] = ICD.count

# Run map
MAP.dir =  "source/FUN_NLP_PheWAS_v2.R"
source(MAP.dir)

dat.ICDcount = list(ID = dat.MAP$patientNum,
                    mat = data.frame(colorectal_cancer = dat.MAP$ICD))
dat.CUIcount = list(ID = dat.MAP$patientNum,
                    mat = data.frame(colorectal_cancer = dat.MAP$CUI))
dat.notecount = list(ID = dat.MAP$patientNum,
                     mat = data.frame(colorectal_cancer =  dat.MAP$note))

require(data.table)
require(bit64)
require(flexmix)
require(Matrix)
require(stats)
MAPfit = MAP_PheWAS_main(dat.icd = dat.ICDcount, dat.nlp = dat.CUIcount, 
                         dat.note = dat.notecount,
                         nm.phe = "colorectal_cancer",
                         p.icd = 0.001, n.icd = 10, 
                         p.nlp = 0.001, n.nlp = 10,
                         yes.con=FALSE, yes.nlp=FALSE)
dat.MAP$MAP = drop(MAPfit$res.all$colorectal_cancer$scores)

write.csv(dat.MAP,
      file = "result/phenotyping/MAP_CRC_ICDday.csv",row.names=F)

# ROC
dat.MAP = read.csv("result/phenotyping/MAP_CRC_ICDday.csv",
                   stringsAsFactors = F)

library(pROC)

# All patients
roc.MAP = roc(dat.MAP$label,dat.MAP$MAP, direction = "<")
roc.ICD = roc(dat.MAP$label,dat.MAP$ICD, direction = "<")
roc.CUI = roc(dat.MAP$label,dat.MAP$CUI, direction = "<")

png("figures/roc_MAP_CRC_ICDday_all.png")
plot(roc.MAP, col = "red")
lines(roc.ICD, col = "blue", lty = 2)
lines(roc.CUI, col = "green", lty = 2)
legend("bottomright", c(paste("MAP, AUC =", format(roc.MAP$auc,digits=0,nsmall=3)),
                        paste("ICD, AUC =", format(roc.ICD$auc,digits=0,nsmall=3)),
                        paste("NLP, AUC =", format(roc.CUI$auc,digits=0,nsmall=3))),
       col = c("red","blue","green"),
       lty = c(1,2,2))
dev.off()

# Filter positive
dat.MAP.positive = dat.MAP[dat.MAP$ICD>0,]
roc.MAP = roc(dat.MAP.positive$label,dat.MAP.positive$MAP, direction = "<")
roc.ICD = roc(dat.MAP.positive$label,dat.MAP.positive$ICD, direction = "<")
roc.CUI = roc(dat.MAP.positive$label,dat.MAP.positive$CUI, direction = "<")


# PPV 95, 90

spec95 = min(roc.MAP$specificities[roc.MAP$specificities >=.95])
spec90 = min(roc.MAP$specificities[roc.MAP$specificities >=.90])

case.pred = coords(roc.MAP,x=c(spec90, spec95), input="specificity", 
                     ret = c("threshold","specificity","sensitivity", "ppv"),
                   transpose = F)
case.pred$total = apply( outer(dat.MAP$MAP, case.pred$threshold, ">="),2,sum)
case.pred$pred = case.pred$total*case.pred$ppv

dat.MAP$MAP.pred.spec90 = (dat.MAP$MAP >= case.pred$threshold[1])
dat.MAP$MAP.pred.spec95 = (dat.MAP$MAP >= case.pred$threshold[2])

write.csv(dat.MAP,
          file = "result/phenotyping/MAP_CRC_ICDday.csv",row.names=F)
