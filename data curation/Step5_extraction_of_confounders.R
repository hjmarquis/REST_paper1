dat =read.csv("data/paper surg/CRC_surg_elig_features.csv",
              stringsAsFactors = F)
dat$first_CRC_dx = as.Date(dat$first_CRC_dx)
dat$trt_date = as.Date(dat$trt_date)

local.dir = ""
dx.dir = paste0(local.dir, 
                "Merck/EHR/Merck_CRC_Diagnosis.csv")
rpdr.dir = paste0(local.dir,
                  "PHS data/Codified dictionary/rpdr_code_codebook_Marquis_i2b2.tsv")

rpdr.map = read.delim(rpdr.dir, stringsAsFactors = F)


# Load the data
ICD.dat = read.delim(dx.dir, sep = '|',
                     stringsAsFactors = F, fileEncoding = "UTF-8-BOM")
colnames(ICD.dat) = c("patientNum", "Date", "Code")
ICD.dat = ICD.dat[ICD.dat$patientNum %in% dat$patientNum,]
ICD.dat$Date = as.Date(ICD.dat$Date)
ICD.dat = ICD.dat[order(ICD.dat$patientNum,ICD.dat$Date),]
ICD.dat$gap = as.numeric(dat$trt_date[match(ICD.dat$patientNum,dat$patientNum)] - ICD.dat$Date, 
                         units = "days")

# variables within 5 years
#------------------------------------------------
ICD.dat = ICD.dat[ICD.dat$gap>=0,]
ICD.dat = ICD.dat[ICD.dat$gap<= 356.2425*5,]

cancer.i2b2 = rpdr.map$i2b2[grep("cancer|neoplasm", rpdr.map$feature_desc,
                                 ignore.case = TRUE)]
except1.i2b2 = intersect(rpdr.map$i2b2[grep("Skin cancer", rpdr.map$feature_desc)], 
                         rpdr.map$i2b2[grep("Squamous cell carcinoma", rpdr.map$feature_desc)])
except2.i2b2 = intersect(rpdr.map$i2b2[grep("Skin cancer", rpdr.map$feature_desc)], 
                         rpdr.map$i2b2[grep("Basal cell carcinoma", rpdr.map$feature_desc)])
except3.i2b2 = rpdr.map$i2b2[rpdr.map$feature_id == "PheCode:180.3"]
CRC.i2b2 = rpdr.map$i2b2[grep("PheCode:153",rpdr.map$feature_id)]


ICD.cancer = ICD.dat[ICD.dat$Code %in% cancer.i2b2,]
ICD.cancer = ICD.cancer[ !ICD.cancer$Code %in% c(except1.i2b2, except2.i2b2,except3.i2b2),]

ICD.cancer.bCRC = ICD.cancer[as.numeric(ICD.cancer$Date - 
                                          dat$first_CRC_dx[match(ICD.cancer$patientNum,
                                                                 dat$patientNum)],
                                        units = "days") <= - 90,]

ICD.cancer.bCRC.sum = data.frame(Code=unique(ICD.cancer.bCRC[,c("Code")]))

ICD.cancer.bCRC.sum$npatient = ICD.cancer.bCRC.sum$ndates = ICD.cancer.bCRC.sum$ncodes =  ICD.cancer.bCRC.sum$med.gap = 00

for(i in 1:nrow(ICD.cancer.bCRC.sum))
{
  tmp.ICD = ICD.cancer.bCRC[ICD.cancer.bCRC$Code==ICD.cancer.bCRC.sum$Code[i],]
  ICD.cancer.bCRC.sum$npatient[i] = length(unique(tmp.ICD$patientNum))
  ICD.cancer.bCRC.sum$ndates[i] = nrow(unique(tmp.ICD[,c("patientNum","Date")]))
  ICD.cancer.bCRC.sum$ncodes[i] = nrow(tmp.ICD)
  ICD.cancer.bCRC.sum$med.gap[i] = median(tmp.ICD$gap)
}
dat$prior_cancer_new = dat$patientNum %in% ICD.cancer.bCRC$patientNum
table(dat$prior_cancer, dat$lap.CPT.filtered, useNA = "ifany")
table(dat$prior_cancer_new, dat$lap.CPT.filtered, useNA = "ifany")

# Update Stage 3 & 4
ICD.stage.aCRC = ICD.cancer[as.numeric(ICD.cancer$Date - 
                                         as.Date(dat$first_CRC_dx[match(ICD.cancer$patientNum,
                                                                        dat$patientNum)]),
                                       units = "days") >= - 90,]
ICD.stage.aCRC = ICD.stage.aCRC[grep("secondary",ICD.stage.aCRC$Description,
                                     ignore.case = T),]
ICD.stage.aCRC = ICD.stage.aCRC[ICD.stage.aCRC$gap>5,]
ICD.stage.aCRC$lymph = grepl("lymph",ICD.stage.aCRC$Description,
                             ignore.case = T)
dat$stage.wICD = dat$stage
tmp.pos = dat$patientNum %in% ICD.stage.aCRC$patientNum[!(ICD.stage.aCRC$lymph)]
dat$stage.wICD[tmp.pos] = pmax(dat$stage[tmp.pos],4,na.rm = T)
tmp.pos = dat$patientNum %in% ICD.stage.aCRC$patientNum[ICD.stage.aCRC$lymph]
dat$stage.wICD[tmp.pos] = pmax(dat$stage[tmp.pos],3,na.rm = T)


# Prohibitive adhesions
i2b2.var = rpdr.map$i2b2[rpdr.map$feature_id == "PheCode:568.1"]

id.var = unique(ICD.dat$patientNum[ICD.dat$Code %in% i2b2.var])
dat$adhesion_new = dat$patientNum %in% id.var
table(dat$adhesion, dat$lap.CPT.filtered)
table(dat$adhesion_new, dat$lap.CPT.filtered)

# Crohn's disease
i2b2.var = rpdr.map$i2b2[rpdr.map$feature_id == "PheCode:555.1"]

id.var = unique(ICD.dat$patientNum[ICD.dat$Code %in% i2b2.var])
dat$crohn_new = dat$patientNum %in% id.var
table(dat$crohn, dat$lap.CPT.filtered)
table(dat$crohn_new, dat$lap.CPT.filtered)

# Familial polyposis: no exact matching (also rare)

# Chronic ulcerative colitis
i2b2.var = rpdr.map$i2b2[rpdr.map$feature_id == "PheCode:555.2"]

id.var = unique(ICD.dat$patientNum[ICD.dat$Code %in% i2b2.var])
dat$cuc_new = dat$patientNum %in% id.var
table(dat$cuc, dat$lap.CPT.filtered)
table(dat$cuc_new, dat$lap.CPT.filtered)

# variables within 1 year
#------------------------------------------------

ICD.dat = ICD.dat[ICD.dat$gap<= 356.2425,]

# Obesity
i2b2.var = rpdr.map$i2b2[grep("PheCode:278",rpdr.map$feature_id) ]

id.var = unique(ICD.dat$patientNum[ICD.dat$Code %in% i2b2.var])
dat$obesity_new = dat$patientNum %in% id.var
table(dat$obesity, dat$lap.CPT.filtered)
table(dat$obesity_new, dat$lap.CPT.filtered)


# variables within 30 days
#------------------------------------------------

ICD.dat = ICD.dat[ICD.dat$gap<= 30,]

# Acutely obstructed colon
i2b2.var = rpdr.map$i2b2[grep("intestinal obstruction",rpdr.map$feature_desc,
                              ignore.case = TRUE) ]
id.var = unique(ICD.dat$patientNum[ICD.dat$Code %in% i2b2.var])
dat$obstruct_new = dat$patientNum %in% id.var
table(dat$obstruct, dat$lap.CPT.filtered)
table(dat$obstruct_new, dat$lap.CPT.filtered)


# Acutely perforated colon
ICD9.var=56983
ICD10.var="ICD10:K63.1"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$perforate = dat$patientNum %in% id.var
table(dat$perforate, dat$lap.CPT.filtered)

write.csv(dat, "data/paper surg/CRC_surg_features.csv",
          row.names = F)