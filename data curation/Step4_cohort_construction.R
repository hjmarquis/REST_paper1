#
# Get the CPT colectomy dates
#=============================================

dat = read.csv("result/CRC_MAP_dem_fu.csv", 
               stringsAsFactors = F)
dat = dat[dat$MAP.pred.spec90,]

dat.CPT = read.csv("data/EHR/CPT_treatment_0715.csv",
                   stringsAsFactors = F)
dat.CPT = dat.CPT[is.element(dat.CPT$PatientNum, dat$patientNum),1:3]
# dat.CPT$Code_num = as.integer(gsub('C','',dat.CPT$Code))
dat.CPT$date = as.Date(sapply(strsplit(dat.CPT$Start_date,' '), '[',1))


# Apply the CPT filter
surg.dict = read.csv("data/filter/old/CPT_colectomy_dictionary_20210201.csv", 
                     stringsAsFactors = F, fileEncoding = "UTF-8-BOM")

surg.dict$LAP = grepl("Laparoscopy", surg.dict$Code_Name)
lap.CPT = surg.dict$Code[surg.dict$LAP]
open.CPT = surg.dict$Code[! surg.dict$LAP]
# surg.CPT = paste('C', c(44140:44160,44204:44212, 45110:45123, 45395:45395), sep='')
surg.CPT = c(lap.CPT, open.CPT)

dat.CPT$surg = is.element(dat.CPT$Code, surg.CPT)
dat.CPT$open = is.element(dat.CPT$Code, open.CPT)
dat.CPT$lap = is.element(dat.CPT$Code, lap.CPT)

dat.CPT = dat.CPT[dat.CPT$surg,]


# Filtered cohorts
first.CPT = dat.CPT[!duplicated(dat.CPT$PatientNum),]
rest.CPT = dat.CPT[dat.CPT$Code!=first.CPT$Code[match(dat.CPT$PatientNum, first.CPT$PatientNum)], ]
rest.CPT = rest.CPT[!duplicated(rest.CPT$PatientNum),]

first.CPT = first.CPT[-which(as.numeric(
  rest.CPT$date[match(first.CPT$PatientNum,rest.CPT$PatientNum)]-
    first.CPT$date, units = "days") == 0),]

match.open = match(dat$patientNum, first.CPT$PatientNum[first.CPT$open])
dat$open.CPT.filtered = !is.na(match.open)
match.lap = match(dat$patientNum, first.CPT$PatientNum[first.CPT$lap])
dat$lap.CPT.filtered = !is.na(match.lap)


dat$open.filtered.date = (first.CPT$date[first.CPT$open])[match.open]
dat$lap.filtered.date = (first.CPT$date[first.CPT$lap])[match.lap]


dat = dat[(dat$open.CPT.filtered | dat$lap.CPT.filtered) & dat$MAP.pred.spec95 , 
          c("patientNum", "MAP", "first_CRC_dx", 
            "birth_date", "death_date","last.date", "first.date", "gender", "race", 
            "open.CPT.filtered","lap.CPT.filtered","open.filtered.date","lap.filtered.date")]

dat$trt_date = dat$open.filtered.date
dat$trt_date[dat$lap.CPT.filtered] = dat$lap.filtered.date[dat$lap.CPT.filtered]

write.csv(dat, file = "data/paper surg/CRC_surg_20210201.csv", 
               row.names = F)

#
# Get stage information
#=============================================

rm(list=objects())

histology.dir = "../REST with Merck/NLP_related/Stage_Histology/"
files.histology =paste(histology.dir,list.files(histology.dir),sep='')
files.histology = files.histology[grep(".txt",files.histology)]

min.gap = -30
max.gap = 365.2425

dat.trt = read.csv("data/paper surg/CRC_surg_20210201.csv", 
                   stringsAsFactors = F)
grp.var = c("open.CPT.filtered","lap.CPT.filtered")
date.var = c("open.filtered.date","lap.filtered.date")
for(var.name in date.var)
{
  dat.trt[,var.name] =as.Date(dat.trt[,var.name])
}

grplist = c("openFlt", "lapFlt")
featurelist = c("T","N","M","stage","hist","meta")

ngrp = length(grplist)
nfeature = length(featurelist)

extract.var = c( outer(grplist, featurelist, 
                       paste,sep='.'))
dat.trt[,extract.var] = NA

for(batch.file in files.histology)
{
  batch.dat = read.delim(batch.file, sep="|",
                         stringsAsFactors = F)
  batch.dat = batch.dat[batch.dat$cuis!="",]
  batch.dat = batch.dat[is.element(batch.dat$patient_num,dat.trt$patientNum),]
  batch.dat = batch.dat[order(batch.dat$patient_num,batch.dat$date),]
  match.pos = match(batch.dat$patient_num,dat.trt$patientNum)
  batch.dat$date = as.Date(as.character(batch.dat$date), "%Y%m%d")
  
  old.id = -1
  for(ipos in 1:nrow(batch.dat))
  {
    datpos = match.pos[ipos]
    if(batch.dat$patient_num[ipos] != old.id)
    {
      old.id = batch.dat$patient_num[ipos]
      conf = rep(0,ngrp*nfeature)
    }
    
    
    extract = strsplit(unlist(strsplit(batch.dat$cuis[ipos], "\\[|\\]|,")),
                       "-|:")
    extract = extract[sapply(extract, length) == 3]
    extract = do.call(rbind, extract)
    
    # Rule: Mention of disease in the month
    extract = extract[extract[,1]!="L", , drop = FALSE]

    
    if(is.null(nrow(extract)))
      next
    
    if(nrow(extract)==0 )
      next
    
    dup.term = duplicated(extract[,1:2])
    if(any(dup.term))
    {
      dup.extract = extract[dup.term,,drop = FALSE]
      extract = extract[!dup.term,,drop = FALSE]
      
      dup.label = apply(dup.extract[,1:2,drop = FALSE],1,paste, collapse = '-')
      extract.label = apply(extract[,1:2,drop = FALSE],1,paste, collapse = '-')
      dup.match = match(dup.label,extract.label)
      
      for(iterm in 1:nrow(dup.extract))
      {
        extract[dup.match[iterm],3] = paste(extract[dup.match[iterm],3],
                                            dup.extract[iterm,3],
                                            sep = ";")
      }
    }
    
    for (iterm in 1:nrow(extract))
    {
      value = unlist(strsplit(extract[iterm,3],";"))
      
      if(extract[iterm,2]=="Clinicalstage")
      {
        extract[iterm,3] = max(as.numeric(gsub("[a-z]", "",gsub("[a-z][0-9]", "", value))))
      }else if(extract[iterm,2]=="TNMstage")
      {
        tval = nval = mval = -1
        for(tmp.val in value)
        {
          tpos = unlist(gregexpr('t',tmp.val))
          tmp.t = as.numeric(substr(tmp.val,tpos+1,tpos+1))
          tval = max(tval,tmp.t, na.rm = T)
          
          npos = unlist(gregexpr('n',tmp.val))
          tmp.n = as.numeric(substr(tmp.val,npos+1,npos+1))
          nval = max(nval,tmp.n, na.rm = T)
          
          mpos = unlist(gregexpr('m',tmp.val))
          tmp.m = as.numeric(substr(tmp.val,mpos+1,mpos+1))
          mval = max(mval,tmp.m, na.rm = T)
        }
        extract[iterm,3] = paste(tval,nval,mval)
      }else{
        extract[iterm,3] = names(sort(table(value)))[1]
      }
    }
    extract[,1] = gsub("M",1,gsub("H",2,extract[,1]))
    
    
    for (igrp in 1:ngrp)
    {
      if(!dat.trt[datpos,grp.var[igrp]])
        next
      
      date.gap = as.numeric(dat.trt[datpos,date.var[igrp]]-batch.dat$date[ipos], units = "days")
      if(date.gap > max.gap | date.gap < min.gap)
        next
      
      # if(date.gap < 0 & !is.na(dat.trt[datpos,ecog.var[igrp]]))
      #   next
      
      for(iterm in 1:nrow(extract))
      {
        conf.diff = as.numeric(extract[iterm,1]) - conf + (date.gap >= 0)
        if(extract[iterm,2]=="Clinicalstage")
        {
          if(conf.diff[3*ngrp+igrp]>0)
          {
            conf[3*ngrp+igrp] = as.numeric(extract[iterm,1])
            dat.trt[datpos, extract.var[3*ngrp+igrp]] = extract[iterm,3]
          }
          next
        }
        
        if(extract[iterm,2]=="Histology")
        {
          if(conf.diff[4*ngrp+igrp]>0)
          {
            conf[4*ngrp+igrp] = as.numeric(extract[iterm,1])
            dat.trt[datpos, extract.var[4*ngrp+igrp]] = extract[iterm,3]
          }
          next
        }
        
        if(extract[iterm,2]=="metastasis")
        {
          if(conf.diff[5*ngrp+igrp]>0)
          {
            conf[5*ngrp+igrp] = as.numeric(extract[iterm,1])
            dat.trt[datpos, extract.var[5*ngrp+igrp]] = extract[iterm,3]
          }
          next
        }
        
        
        tnm = unlist(strsplit(extract[iterm,3], " "))
        
        if(tnm[1] >=0 &
           conf.diff[igrp]>0)
        {
          conf[igrp] = as.numeric(extract[iterm,1])
          dat.trt[datpos, extract.var[igrp]] = tnm[1]
        }
        
        if(tnm[2] >=0 &
           conf.diff[ngrp+igrp]>0)
        {
          conf[ngrp+igrp] = as.numeric(extract[iterm,1])
          dat.trt[datpos, extract.var[ngrp+igrp]] = tnm[2]
        }
        
        if(tnm[3] >=0 &
           conf.diff[2*ngrp+igrp]>0)
        {
          conf[2*ngrp+igrp] = as.numeric(extract[iterm,1])
          dat.trt[datpos, extract.var[2*ngrp+igrp]] = tnm[3]
        }
      }
    }
  }
}

write.csv(dat.trt, "data/paper surg/CRC_surg_20210201_stage.csv",
          row.names = F)

#
# Get BMI information
#=============================================

rm(list=objects())

# Trim the data
dat.trt = read.csv("data/paper surg/CRC_surg_20210201_stage.csv", stringsAsFactors = F)


dat.trt$trt_date = dat.trt$open.filtered.date
dat.trt$trt_date[dat.trt$lap.CPT.filtered] = dat.trt$lap.filtered.date[dat.trt$lap.CPT.filtered]

write.csv(dat.trt[,c("patientNum","trt_date")],
          file = "data/paper surg/CRC_surg_20210201_patientNum_date.csv")

dat.trt$stage = dat.trt$openFlt.stage
dat.trt$stage[dat.trt$lap.CPT.filtered] = dat.trt$lapFlt.stage[dat.trt$lap.CPT.filtered]

stageM = dat.trt$openFlt.M
stageM[dat.trt$lap.CPT.filtered] = dat.trt$lapFlt.M[dat.trt$lap.CPT.filtered]
dat.trt$stage[which(stageM>=0.5)] = 4

stageN = dat.trt$openFlt.N
stageN[dat.trt$lap.CPT.filtered] = dat.trt$lapFlt.N[dat.trt$lap.CPT.filtered]
dat.trt$stage[which(stageN>=0.5)] = pmax(3,dat.trt$stage[which(stageN>=0.5)], na.rm = T)

stageT = dat.trt$openFlt.T
stageT[dat.trt$lap.CPT.filtered] = dat.trt$lapFlt.T[dat.trt$lap.CPT.filtered]
dat.trt$stage[which(stageT>=2.5)] = pmax(2,dat.trt$stage[which(stageT>=2.5)], na.rm = T)
dat.trt$stage[which(stageT>=0.5)] = pmax(1,dat.trt$stage[which(stageT>=0.5)], na.rm = T)
dat.trt$stage[which(stageT<0.5)] = pmax(0,dat.trt$stage[which(stageT<0.5)], na.rm = T)

dat.trt = dat.trt[,c("patientNum", "birth_date", "death_date", 
                     "first.date","last.date", "first_CRC_dx", 
                     "gender", "race", 
                     "lap.CPT.filtered","trt_date","stage")]

# Get the BMI from codified data
#---------------------------------------------------------------------------------
dat.bmi = read.delim("NLP_related/BMI_data/BMI_colorectal.csv",
                     sep = '|', stringsAsFactors = F, fileEncoding = "UTF-8-BOM")
dat.bmi = dat.bmi[dat.bmi$PatientNum %in% dat.trt$patientNum, ]
dat.bmi$Start_date = as.Date(dat.bmi$Start_date)
dat.trt$trt_date = as.Date(dat.trt$trt_date)

has.bmi = unique(dat.bmi$PatientNum[dat.bmi$Code_Name=="BMI"])
has.wgt.hgt = intersect(unique(dat.bmi$PatientNum[dat.bmi$Code_Name=="Weight"]),
                        unique(dat.bmi$PatientNum[dat.bmi$Code_Name=="Height"]))
has.either = unique(c(has.bmi,has.wgt.hgt))

dat.bmi = dat.bmi[order(dat.bmi$PatientNum, dat.bmi$Start_date),]
dat.bmi$gap = as.numeric(dat.trt$trt_date[match(dat.bmi$PatientNum, dat.trt$patientNum)] - dat.bmi$Start_date,
                         units = "days")
# Clean the BMI data
dat.wgt = dat.bmi[dat.bmi$Code_Name=="Weight",]
dat.hgt = dat.bmi[dat.bmi$Code_Name=="Height",]
dat.bmi = dat.bmi[dat.bmi$Code_Name=="BMI",]

inch.to.cm = 2.54
dat.hgt = dat.hgt[(as.numeric(dat.hgt$Start_date-
                                as.Date(dat.trt$birth_date[match(dat.hgt$PatientNum,dat.trt$patientNum)]),
                              units = "days")/356.2425) >= 18,]
dat.hgt$hgt.clean = as.numeric(dat.hgt$Result)
dat.hgt = dat.hgt[!is.na(dat.hgt$hgt.clean),]

tmp.hgt = dat.hgt$hgt.clean * inch.to.cm
inch.div.254 = which((tmp.hgt< 100) & (tmp.hgt>=30) &
                       ((abs(tmp.hgt-round(tmp.hgt,2))<0.001)|
                          (abs(tmp.hgt-round(tmp.hgt,1))<0.01)|
                          (abs(tmp.hgt-round(tmp.hgt))<0.1)
                       ))
ft.div.254 = which((tmp.hgt< 8) & 
                     (abs(tmp.hgt-round(tmp.hgt,1))<0.01))
dat.hgt$hgt.clean[c(inch.div.254,ft.div.254)] = tmp.hgt[c(inch.div.254,ft.div.254)]

ft.inch = which((dat.hgt$hgt.clean <8) &
                  (dat.hgt$hgt.clean > 4))
dat.hgt$hgt.clean[ft.inch] = dat.hgt$hgt.clean[ft.inch]*10 + round(dat.hgt$hgt.clean[ft.inch])*2
dat.hgt = dat.hgt[dat.hgt$hgt.clean>0,]

pat.first = which(!duplicated(dat.hgt$PatientNum))
pat.last = c(pat.first[-1]-1, nrow(dat.hgt))
pat.rep = pat.last-pat.first +1
pat.first = rep(pat.first,pat.rep)
pat.last = rep(pat.last,pat.rep)

dat.hgt$hgt.clean2 = dat.hgt$hgt.clean
for (i in 1:nrow(dat.hgt))
{
  if(dat.hgt$hgt.clean[i]>100)
    next
  tmp.hgt = dat.hgt$hgt.clean[i] + round(dat.hgt$hgt.clean[i]/10)*2
  other.pos = (pat.first[i]:pat.last[i])[dat.hgt$Code[pat.first[i]:pat.last[i]]!=4994]
  if(any((abs(dat.hgt$hgt.clean[other.pos]-tmp.hgt)/tmp.hgt)<0.05))
    dat.hgt$hgt.clean2[i] = tmp.hgt
  if(any((abs(dat.hgt$hgt.clean[other.pos]/inch.to.cm-tmp.hgt)/tmp.hgt)<0.05))
    dat.hgt$hgt.clean2[i] = tmp.hgt
}

hgt.inch = dat.hgt$hgt.clean < 100
dat.hgt$hgt.clean[hgt.inch] = dat.hgt$hgt.clean[hgt.inch] * inch.to.cm
hgt.inch = dat.hgt$hgt.clean2 < 100
dat.hgt$hgt.clean2[hgt.inch] = dat.hgt$hgt.clean2[hgt.inch] * inch.to.cm

yo = unique(dat.hgt$PatientNum[dat.hgt$hgt.clean2<135])
yo = unique(dat.hgt$PatientNum[dat.hgt$hgt.clean2>220])
i=0
i=i+1;dat.hgt[dat.hgt$PatientNum==yo[i],]
summary(dat.hgt$hgt.clean)
summary(dat.hgt$hgt.clean2)
dat.hgt[dat.hgt$hgt.clean2<=230,]

dat.trt$hgt = NA
map.pos = match(dat.hgt$PatientNum,dat.trt$patientNum)
for(i in unique(pat.first))
{
  dat.trt$hgt[map.pos[i]] = median(dat.hgt$hgt.clean2[pat.first[i]:pat.last[i]])
}

dat.wgt$Result = gsub(" s\\.| s", '',dat.wgt$Result)
dat.wgt$wgt.clean = as.numeric(dat.wgt$Result)
dat.wgt = dat.wgt[!is.na(dat.wgt$wgt.clean),]
dat.wgt = dat.wgt[dat.wgt$wgt.clean >0,]

mis.dec = dat.wgt$wgt.clean >= 1000
dat.wgt$wgt.clean[mis.dec] = dat.wgt$wgt.clean[mis.dec]/10


lb.to.kg = 1/2.205
wgt.lb = dat.wgt$Code!=4993
dat.wgt$wgt.clean[wgt.lb] = dat.wgt$wgt.clean[wgt.lb]*lb.to.kg

med.win.max = 365.2425
med.win.min = -365.2425

i.id = -1
dat.wgt$wgt.med365 = NA
for (i in 1:nrow(dat.wgt))
{
  if(dat.wgt$PatientNum[i]!=i.id)
  {
    i.id = dat.wgt$PatientNum[i]
    i.start = i.end = i
  }
  
  while ((as.numeric(dat.wgt$Start_date[i]-dat.wgt$Start_date[i.start],
                     units = "days")>med.win.max))
  {
    i.start = i.start +1
  }
  
  next.end = min(nrow(dat.wgt),i.end+1)
  while ((dat.wgt$PatientNum[next.end]==i.id)&
         (as.numeric(dat.wgt$Start_date[i]-dat.wgt$Start_date[next.end],
                     units = "days")>med.win.min))
  {
    i.end = next.end
    next.end = next.end+1
    if(next.end>nrow(dat.wgt))
      break
  }
  dat.wgt$wgt.med365[i] = median(dat.wgt$wgt.clean[i.start:i.end])
}

summary(dat.wgt$wgt.clean[dat.wgt$Code!=4993])
yo = unique(dat.wgt$PatientNum[dat.wgt$wgt.med365<20])
i=0
i=i+1;dat.wgt[dat.wgt$PatientNum==yo[i],]

summary(dat.wgt$wgt.med365)

dat.wgt$bmi.med365 = dat.wgt$wgt.med365/((dat.trt$hgt[match(dat.wgt$PatientNum,dat.trt$patientNum)]/100)^2)

summary(dat.wgt$bmi.med365)

dat.bmi$bmi.clean = as.numeric(dat.bmi$Result)

yo = unique(dat.bmi$PatientNum[dat.bmi$bmi.clean<10])
yo = unique(dat.bmi$PatientNum[dat.bmi$bmi.clean>70])
i=0
i=i+1;dat.bmi[dat.bmi$PatientNum==yo[i],]

med.win.max = 365.2425
med.win.min = -365.2425

i.id = -1
dat.bmi$bmi.med365 = NA
for (i in 1:nrow(dat.bmi))
{
  if(dat.bmi$PatientNum[i]!=i.id)
  {
    i.id = dat.bmi$PatientNum[i]
    i.start = i.end = i
  }
  
  while ((as.numeric(dat.bmi$Start_date[i]-dat.bmi$Start_date[i.start],
                     units = "days")>med.win.max))
  {
    i.start = i.start +1
  }
  
  next.end = min(nrow(dat.bmi),i.end+1)
  while ((dat.bmi$PatientNum[next.end]==i.id)&
         (as.numeric(dat.bmi$Start_date[i]-dat.bmi$Start_date[next.end],
                     units = "days")>med.win.min))
  {
    i.end = next.end
    next.end = next.end+1
    if(next.end>nrow(dat.bmi))
      break
  }
  dat.bmi$bmi.med365[i] = median(dat.bmi$bmi.clean[i.start:i.end])
}

summary(dat.bmi$bmi.med365)

# Get BMI
dat.bmi.clean = rbind(dat.bmi[,c("PatientNum", "Start_date","bmi.med365", "gap")],
                      dat.wgt[,c("PatientNum", "Start_date","bmi.med365", "gap")])
dat.bmi.clean = dat.bmi.clean[!is.na(dat.bmi.clean$bmi.med365),]
dat.bmi.clean = dat.bmi.clean[order(dat.bmi.clean$PatientNum, 
                                    dat.bmi.clean$bmi.med365<15,
                                    ((dat.bmi.clean$gap<0) | (dat.bmi.clean$gap>365.2425)),
                                    ((dat.bmi.clean$gap<= -30) | (dat.bmi.clean$gap>365.2425)),
                                    (dat.bmi.clean$gap<0),
                                    abs(dat.bmi.clean$gap)),]
dat.trt$bmi.direct = NA
dat.trt$bmi.direct.gap = NA

map.pos = match(dat.bmi.clean$PatientNum, dat.trt$patientNum)
best.pos = !duplicated(dat.bmi.clean$PatientNum)
dat.trt$bmi.direct[map.pos[best.pos]] = dat.bmi.clean$bmi.med365[best.pos]
dat.trt$bmi.direct.gap[map.pos[best.pos]] = dat.bmi.clean$gap[best.pos]


summary(dat.trt$bmi.direct)
summary(dat.trt$bmi.direct.gap)

write.csv(dat.bmi.clean, file="data/paper surg/RCTsurg_BMI_clean_Structure.csv",
          row.names = F)

write.csv(dat.trt, file = "data/paper surg/CRC_surg_20210201_stage_bmi.csv",
          row.names = F)

dat.trt$nobmi = is.na(dat.trt$bmi.direct)
dat.trt$nobmi_at_surg = F
dat.trt$nobmi_at_surg[which(dat.trt$bmi.direct.gap<= -30)] = T
dat.trt$nobmi_at_surg[which(dat.trt$bmi.direct.gap> 365.2425)] = T

dat.missbmi = dat.trt[dat.trt$nobmi | dat.trt$nobmi_at_surg,
                      c("patientNum","trt_date",
                        "bmi.direct", "bmi.direct.gap",
                        "nobmi", "nobmi_at_surg")]

write.csv(dat.missbmi, file = "data/paper surg/CRC_surg_20210201_misbmi.csv",
          row.names = F)


# NLP 
#-----------------------------------------------

rm(list = objects())

nlp.read = strsplit(read.delim("../REST with Merck/NLP_related/BMI_data/EXTEND_BMI.csv",
                               stringsAsFactors = F,
                               header = F)$V1, "\\|")
nlp.read = nlp.read[sapply(nlp.read, length)>3]

nextract = (sapply(nlp.read, length)-3)/2

dat.bmi = data.frame(PatientNum = rep(NA,sum(nextract)),
                     Start_date = as.Date(NA),
                     Code_Name = "",
                     Result = NA,
                     Reference_Units = "",
                     stringsAsFactors = F)

map.pos = rep(1:length(nlp.read),nextract)
extract.pos = unlist(sapply(sapply(nextract, ":", 1),rev))

for (i in 1:nrow(dat.bmi))
{
  dat.bmi$PatientNum[i] = nlp.read[[map.pos[i]]][1]
  dat.bmi$Start_date[i] = as.Date(nlp.read[[map.pos[i]]][3])
  dat.bmi$Code_Name[i] = nlp.read[[map.pos[i]]][2+2*extract.pos[i]]
  dat.bmi$Result[i] = gsub(" |[a-z]|[A-Z]", '',nlp.read[[map.pos[i]]][3+2*extract.pos[i]])
  dat.bmi$Reference_Units[i] = gsub(" |[0-9]|\\.",'',nlp.read[[map.pos[i]]][3+2*extract.pos[i]])
}

dat.bmi = dat.bmi[order(dat.bmi$PatientNum, dat.bmi$Start_date),]

write.csv(dat.bmi, file = "NLP_related/BMI_data/EXTEND_BMI_cleaned.csv",
          row.names = F)

dat.trt = read.csv("data/paper surg/CRC_surg_20210201_stage_bmi.csv",
                   stringsAsFactors = F)
dat.trt$trt_date = as.Date(dat.trt$trt_date)
dat.bmi$gap = as.numeric(dat.trt$trt_date[match(dat.bmi$PatientNum, dat.trt$patientNum)] - dat.bmi$Start_date,
                         units = "days")

# Standardize the units
#-------------------------------

# Height
dat.hgt = dat.bmi[dat.bmi$Code_Name=="H",]
dat.hgt = dat.hgt[!(dat.hgt$Reference_Units %in% c("f","m")),]

dat.hgt$hgt.clean = NA

tmp.pos = dat.hgt$Reference_Units=="cm"
dat.hgt$hgt.clean[tmp.pos] = as.numeric(dat.hgt$Result[tmp.pos])
tmp.pos = grep("^in", dat.hgt$Reference_Units)
dat.hgt$hgt.clean[tmp.pos] = as.numeric(substr(dat.hgt$Result[tmp.pos],1,2))*2.54
tmp.pos =  grep("^f", dat.hgt$Reference_Units)
dat.hgt$hgt.clean[tmp.pos] = as.numeric(substr(dat.hgt$Result[tmp.pos],1,1))*12*2.54
tmp.pos =  (grepl("^f", dat.hgt$Reference_Units) & grepl("in", dat.hgt$Reference_Units) &
              (nchar(dat.hgt$Result)>=2))
dat.hgt$hgt.clean[tmp.pos] = (dat.hgt$hgt.clean[tmp.pos] + 
                                as.numeric(substr(dat.hgt$Result[tmp.pos],2,2))*2.54)
dat.hgt = dat.hgt[dat.hgt$hgt.clean>=140,]

pat.first = which(!duplicated(dat.hgt$PatientNum))
pat.last = c(pat.first[-1]-1, nrow(dat.hgt))
pat.rep = pat.last-pat.first +1
pat.first = rep(pat.first,pat.rep)
pat.last = rep(pat.last,pat.rep)

dat.trt$hgt.nlp = NA
map.pos = match(dat.hgt$PatientNum,dat.trt$patientNum)
for(i in unique(pat.first))
{
  dat.trt$hgt.nlp[map.pos[i]] = median(dat.hgt$hgt.clean[pat.first[i]:pat.last[i]])
}

dat.trt$hgt.comb = dat.trt$hgt
tmp.pos = is.na(dat.trt$hgt)
dat.trt$hgt.comb[tmp.pos] = dat.trt$hgt.nlp[tmp.pos]
tmp.pos = which(abs(170-dat.trt$hgt) > abs(170-dat.trt$hgt.nlp))
dat.trt$hgt.comb[tmp.pos] = dat.trt$hgt.nlp[tmp.pos]

# Weight
dat.wgt = dat.bmi[dat.bmi$Code_Name=="W",]
dat.wgt = dat.wgt[-grep("^g|to",dat.wgt$Reference_Units),]

dat.wgt$wgt.clean = NA

tmp.pos = grep("^k", dat.wgt$Reference_Units)
dat.wgt$wgt.clean[tmp.pos] = as.numeric(dat.wgt$Result[tmp.pos])
dat.wgt$wgt.clean[-tmp.pos] = as.numeric(dat.wgt$Result[-tmp.pos])/2.205

med.win.max = 365.2425
med.win.min = -365.2425

i.id = -1
dat.wgt$wgt.med365 = NA
for (i in 1:nrow(dat.wgt))
{
  if(dat.wgt$PatientNum[i]!=i.id)
  {
    i.id = dat.wgt$PatientNum[i]
    i.start = i.end = i
  }
  
  while ((as.numeric(dat.wgt$Start_date[i]-dat.wgt$Start_date[i.start],
                     units = "days")>med.win.max))
  {
    i.start = i.start +1
  }
  
  next.end = min(nrow(dat.wgt),i.end+1)
  while ((dat.wgt$PatientNum[next.end]==i.id)&
         (as.numeric(dat.wgt$Start_date[i]-dat.wgt$Start_date[next.end],
                     units = "days")>med.win.min))
  {
    i.end = next.end
    next.end = next.end+1
    if(next.end>nrow(dat.wgt))
      break
  }
  dat.wgt$wgt.med365[i] = median(dat.wgt$wgt.clean[i.start:i.end])
}

dat.wgt$bmi.med365 = dat.wgt$wgt.med365/((dat.trt$hgt.comb[match(dat.wgt$PatientNum,dat.trt$patientNum)]/100)^2)

# BMI
dat.bmi = dat.bmi[dat.bmi$Code_Name=="BMI",]
dat.bmi$bmi.clean = as.numeric(dat.bmi$Result)

med.win.max = 365.2425
med.win.min = -365.2425

i.id = -1
dat.bmi$bmi.med365 = NA
for (i in 1:nrow(dat.bmi))
{
  if(dat.bmi$PatientNum[i]!=i.id)
  {
    i.id = dat.bmi$PatientNum[i]
    i.start = i.end = i
  }
  
  while ((as.numeric(dat.bmi$Start_date[i]-dat.bmi$Start_date[i.start],
                     units = "days")>med.win.max))
  {
    i.start = i.start +1
  }
  
  next.end = min(nrow(dat.bmi),i.end+1)
  while ((dat.bmi$PatientNum[next.end]==i.id)&
         (as.numeric(dat.bmi$Start_date[i]-dat.bmi$Start_date[next.end],
                     units = "days")>med.win.min))
  {
    i.end = next.end
    next.end = next.end+1
    if(next.end>nrow(dat.bmi))
      break
  }
  dat.bmi$bmi.med365[i] = median(dat.bmi$bmi.clean[i.start:i.end])
}

# Get BMI
dat.bmi.clean = rbind(dat.bmi[,c("PatientNum", "Start_date","bmi.med365", "gap")],
                      dat.wgt[,c("PatientNum", "Start_date","bmi.med365", "gap")])
dat.bmi.clean = dat.bmi.clean[!is.na(dat.bmi.clean$bmi.med365),]
dat.bmi.clean = dat.bmi.clean[order(dat.bmi.clean$PatientNum, 
                                    dat.bmi.clean$bmi.med365<15,
                                    ((dat.bmi.clean$gap<0) | (dat.bmi.clean$gap>365.2425)),
                                    ((dat.bmi.clean$gap<= -30) | (dat.bmi.clean$gap>365.2425)),
                                    (dat.bmi.clean$gap<0),
                                    abs(dat.bmi.clean$gap)),]
dat.trt$bmi.nlp = NA
dat.trt$bmi.nlp.gap = NA

dat.bmi.clean = dat.bmi.clean[dat.bmi.clean$PatientNum %in% dat.trt$patientNum,]
map.pos = match(dat.bmi.clean$PatientNum, dat.trt$patientNum)
best.pos = !duplicated(dat.bmi.clean$PatientNum)
dat.trt$bmi.nlp[map.pos[best.pos]] = dat.bmi.clean$bmi.med365[best.pos]
dat.trt$bmi.nlp.gap[map.pos[best.pos]] = dat.bmi.clean$gap[best.pos]

write.csv(dat.bmi.clean, file="data/paper surg/RCTsurg_BMI_clean_NLP.csv",
          row.names = F)

# Combine with Structured data
dat.bmi.clean = rbind(read.csv("data/paper surg/RCTsurg_BMI_clean_NLP.csv",
                               stringsAsFactors = F),
                      read.csv("data/paper surg/RCTsurg_BMI_clean_Structure.csv",
                               stringsAsFactors = F))
dat.bmi.clean = dat.bmi.clean[order(dat.bmi.clean$PatientNum, 
                                    dat.bmi.clean$bmi.med365<15,
                                    ((dat.bmi.clean$gap<0) | (dat.bmi.clean$gap>365.2425)),
                                    ((dat.bmi.clean$gap<= -30) | (dat.bmi.clean$gap>365.2425)),
                                    (dat.bmi.clean$gap<0),
                                    abs(dat.bmi.clean$gap)),]

dat.trt$bmi.comb = NA
dat.trt$bmi.comb.gap = NA

map.pos = match(dat.bmi.clean$PatientNum, dat.trt$patientNum)
best.pos = !duplicated(dat.bmi.clean$PatientNum)
dat.trt$bmi.comb[map.pos[best.pos]] = dat.bmi.clean$bmi.med365[best.pos]
dat.trt$bmi.comb.gap[map.pos[best.pos]] = dat.bmi.clean$gap[best.pos]

write.csv(dat.trt, file = "data/paper surg/CRC_surg_20210201_stage_bmi.csv",
          row.names = F)

#
#  Other features
#==================================================

rm(list=objects())

dat = read.csv("data/paper surg/CRC_surg_20210201_stage_bmi.csv", 
               stringsAsFactors = F)
dat$trt_date = as.Date(dat$trt_date)
dat$last.date = as.Date(dat$last.date)
dat$first.date = as.Date(dat$first.date)
dat$first_CRC_dx = as.Date(dat$first_CRC_dx)

dat$history = as.numeric(dat$trt_date - dat$first.date, units = "days")/365.2425
# dat$history.dx = as.numeric(dat$first_CRC_dx - dat$first.date, units = "days")/365.2425

#
#  ICD data
#====================================================
local.dir = ""
dx.dir = paste0(local.dir, 
                "Merck/EHR/Merck_CRC_Diagnosis.csv")

# Load the data
ICD.dat = read.delim(dx.dir, sep = '|',
                     stringsAsFactors = F, fileEncoding = "UTF-8-BOM")
colnames(ICD.dat) = c("patientNum", "Date", "Code")
ICD.dat = ICD.dat[ICD.dat$patientNum %in% dat$patientNum,]
ICD.dat$Date = as.Date(ICD.dat$Date)
ICD.dat = ICD.dat[order(ICD.dat$patientNum,ICD.dat$Date),]
ICD.dat$gap = as.numeric(dat$trt_date[match(ICD.dat$patientNum,dat$patientNum)] - ICD.dat$Date, 
                         units = "days")

# Refined disease: search one year before the treatment
sgCRC.ICD10 = c(paste("ICD10:C18",c(1,0,2:7),sep='.'), "ICD10:C19", "ICD10:C20")
sgCRC.ICD9 = c( 1535, 1534, 1536, 1530, 1531, 1537, 1532, 1533, 1540,
                154)

# seg = c("appendix","caecum",
#         "ascending colon", "hepatic flexure",
#         "transverse colon", "splenic flexure",
#         "descending colon", "sigmoid colon", "rectosigmoid junction", 
#         "rectum", 
#         "multiple", "unknown")
seg = c("Right","Right",
        "Right", "Right",
        "transverse colon", "Left",
        "Left", "sigmoid", "sigmoid", 
        "rectum", 
        "multiple", "unknown")
seg.uniq = unique(seg)

ICD.CRCseg = ICD.dat[ grep(paste('^',c(sgCRC.ICD10, sgCRC.ICD9),sep='',collapse = '|'),
                           ICD.dat$Code),]
ICD.CRCseg = ICD.CRCseg[(ICD.CRCseg$gap<30) & (ICD.CRCseg$gap>= -7),]
unique(ICD.CRCseg$Code)

table(dat$patientNum %in% ICD.CRCseg$patientNum, dat$lap.CPT.filtered)

ICD.CRCseg = ICD.CRCseg[order(ICD.CRCseg$patientNum,
                              ICD.CRCseg$gap<0,
                              abs(ICD.CRCseg$gap)),]
best.gap = ICD.CRCseg[!duplicated(ICD.CRCseg$patientNum),]
ICD.CRCseg$best.gap = best.gap$gap[match(ICD.CRCseg$patientNum, best.gap$patientNum)]
ICD.CRCseg = ICD.CRCseg[ICD.CRCseg$gap==ICD.CRCseg$best.gap,]
dat$best.gap = best.gap$gap[match(dat$patientNum, best.gap$patientNum)]

dat$CRC_seg_surgday = factor(NA, levels = seg.uniq)
dat$CRC_seg2_surgday = F
for (i in 1:nrow(dat))
{
  tmp.ICD = ICD.CRCseg[ICD.CRCseg$patientNum==dat$patientNum[i],]
  tmp.ICD$Code = gsub(paste('^',1541:1549,sep='',collapse = '|'),
                      "154",tmp.ICD$Code)
  tmp.seg.pos = which(table(factor(tmp.ICD$Code, c(sgCRC.ICD10, sgCRC.ICD9),
                                   labels = rep(seg[1:10],2)))>0)
  if(any(names(tmp.seg.pos)=="sigmoid"))
  {
    tmp.seg.pos = tmp.seg.pos[setdiff(names(tmp.seg.pos),c("rectum"))]
  }else if(any(names(tmp.seg.pos)=="Left"))
  {
    tmp.seg.pos = tmp.seg.pos[setdiff(names(tmp.seg.pos),c("sigmoid"))]
  }
  tmp.seg = seg.uniq[tmp.seg.pos]
  if(length(tmp.seg.pos)>= 2)
  {
    tmp.seg = "multiple"
  }else if(length(tmp.seg) == 0)
  {
    if(nrow(tmp.ICD)>0)
      stop("check")
    tmp.seg = "unknown"
  }
  
  dat$CRC_seg_surgday[i] = tmp.seg
}

table(dat$CRC_seg_surgday,  dat$lap.CPT.filtered)

# variables within 5 years
#------------------------------------------------
ICD.dat = ICD.dat[ICD.dat$gap>=0,]

ICD.dat = ICD.dat[ICD.dat$gap<= 356.2425*5,]

# Other cancer
cancer.ICD10 = c(paste("ICD10:C",c(0,10:96),sep=''),
                 paste("ICD10:D",c(0,37:49),sep=''))
cancer.ICD9 = c(140:209, 230:239)
except.ICD10 = c("ICD10:C44.92","ICD10:C44.91","ICD10:D06.9")
except.ICD9 = c(17392, 17391, 2331)
CRC.ICD10.nos =  c("ICD10:C18","ICD10:C18.8","ICD10:C18.9")
CRC.ICD9.nos =  c(153,1538,1539)
ICD.cancer = ICD.dat[ grep(paste('^',c(cancer.ICD10, cancer.ICD9),sep='',collapse = '|'),
                           ICD.dat$Code),]
ICD.cancer = ICD.cancer[ -grep(paste('^',c(except.ICD10, except.ICD9),sep='',collapse = '|'),
                               ICD.dat$Code),]
# write.csv(ICD.cancer, file = "data/review/RCTsurg-ICD-cancer.csv",
#           row.names = F)
# dat.trt = read.csv("result/CRC_surg_filtered_v4.csv",
#                    stringsAsFactors = F)

ICD.cancer.bCRC = ICD.cancer[as.numeric(ICD.cancer$Date - 
                                          as.Date(dat$first_CRC_dx[match(ICD.cancer$patientNum,
                                                                         dat$patientNum)]),
                                        units = "days") <= - 90,]
# ICD.cancer.bCRC = ICD.cancer.bCRC[-grep("beni",ICD.cancer.bCRC$Description,
#                                         ignore.case = T),]
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

# write.csv(ICD.cancer.bCRC.sum,file =  "data/review/RCTsurg-ICD-prior-cancer-summary.csv",
#           row.names = F)
dat$prior_cancer = dat$patientNum %in% ICD.cancer.bCRC$patientNum
table(dat$prior_cancer, dat$lap.CPT.filtered, useNA = "ifany")

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
ICD9.var="5680"
ICD10.var="ICD10:K66.0"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$adhesion = dat$patientNum %in% id.var
table(dat$adhesion, dat$lap.CPT.filtered)

# Crohn's disease
ICD9.var="555"
ICD10.var="ICD10:K50"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$crohn = dat$patientNum %in% id.var
table(dat$crohn, dat$lap.CPT.filtered)

# Familial polyposis: no exact matching (also rare)

# Chronic ulcerative colitis
ICD9.var="5566"
ICD10.var="ICD10:K51"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$cuc = dat$patientNum %in% id.var
table(dat$cuc, dat$lap.CPT.filtered)

# variables within 1 year
#------------------------------------------------

ICD.dat = ICD.dat[ICD.dat$gap<= 356.2425,]

# Obesity
ICD9.var=278
ICD10.var="ICD10:E66"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$obesity = dat$patientNum %in% id.var
table(dat$obesity, dat$lap.CPT.filtered)


# variables within 30 days
#------------------------------------------------

ICD.dat = ICD.dat[ICD.dat$gap<= 30,]

# Acutely obstructed colon
ICD9.var=c(56081, 5609) 
ICD10.var=paste("ICD10:K56",c(52,609,699),sep='.')

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$obstruct = dat$patientNum %in% id.var
table(dat$obstruct, dat$lap.CPT.filtered)

# Acutely perforated colon
ICD9.var=56983
ICD10.var="ICD10:K63.1"

id.var = unique(ICD.dat$patientNum[grep(paste('^',c(ICD9.var, ICD10.var),sep='',collapse = '|'),
                                        ICD.dat$Code)])
dat$perforate = dat$patientNum %in% id.var
table(dat$perforate, dat$lap.CPT.filtered)

write.csv(dat, "data/paper surg/CRC_surg_elig_features.csv",
          row.names = F)
