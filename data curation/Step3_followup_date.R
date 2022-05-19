dat = read.csv("result/CRC_MAP_dem.csv", 
               stringsAsFactors = F)
dat = dat[dat$MAP.pred.spec90,]
dat$last.ICD = dat$first.ICD = dat$last.NLP = dat$first.NLP = as.Date(NA)

# Get first date and last date NLP
batch.dir = "EHR/NLP-global"
batch.list = sort(list.files(batch.dir))
nbatch = length(batch.list)


dat.pos = 1
last.date = as.Date(NA)

for(ibatch in 1:nbatch)
{
  nlp.dat = read.delim(paste(batch.dir,batch.list[ibatch], sep='/'), sep = '|',
                       stringsAsFactors = F,
                       header = F)
  nlp.dat$V3 = as.Date(as.character(nlp.dat$V3), "%Y%m%d")
  for(ipos in 1:nrow(nlp.dat))
  {
    old.dat.pos = dat.pos
    while (dat$patientNum[dat.pos] < nlp.dat$V1[ipos]) 
    {
      dat.pos = dat.pos + 1
    }
    if(is.na(dat$first.NLP[dat.pos]))
    {
      dat$first.NLP[dat.pos] = nlp.dat$V3[ipos]
      dat$last.NLP[old.dat.pos] = last.date
    }
    last.date = nlp.dat$V3[ipos]
  }
}
if(dat.pos <= nrow(dat))
{
  dat$last.NLP[dat.pos] = last.date
}else{
  dat$last.NLP[old.dat.pos] = last.date
}

# Get first date and last date ICD
ICD.dat = read.csv("data/Diagnosis_codes.csv",
                   stringsAsFactors = F, header = F, 
                   fileEncoding = "UTF-8-BOM")
colnames(ICD.dat)[1:4] = c("patientNum", "Date", "Code", "Description")
ICD.dat = ICD.dat[,1:4]
ICD.dat$Date = as.Date(ICD.dat$Date)
ICD.dat = ICD.dat[order(ICD.dat$patientNum,ICD.dat$Date),]

dat.pos = 1
last.date = as.Date(NA)

for (ipos in 1:nrow(ICD.dat))
{
  old.dat.pos = dat.pos
  while (dat$patientNum[dat.pos] <ICD.dat$patientNum[ipos]) 
  {
    dat.pos = dat.pos + 1
  }
  if(is.na(dat$first.ICD[dat.pos]))
  {
    dat$first.ICD[dat.pos] = ICD.dat$Date[ipos]
    dat$last.ICD[old.dat.pos] = last.date
  }
  last.date = ICD.dat$Date[ipos]
}
if(dat.pos <= nrow(dat))
{
  dat$last.ICD[dat.pos] = last.date
}else{
  dat$last.ICD[old.dat.pos] = last.date
}

dat$first.date = pmin(dat$first.ICD,dat$first.NLP, na.rm = TRUE)
dat$last.date = pmax(dat$last.ICD,dat$last.NLP, na.rm = TRUE)

write.csv(dat, file = "result/CRC_MAP_dem_fu.csv", 
          row.names = F)
