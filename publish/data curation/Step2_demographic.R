
dat.dem = read.csv("data/EHR/Demographic.csv",
                   stringsAsFactors = F)


source("source/cleaning/MMDDYY_to_date.R")

dat.dem$Date_of_Birth = MMDDYY.to.date(dat.dem$Date_of_Birth, dat.dem$Age)
dat.dem$Date_of_Death = MMDDYY.to.date(dat.dem$Date_of_Death)

# write into grand summary
#--------------------------------------------------------
dat = read.csv("result/phenotyping/MAP_CRC_ICDday.csv",
               stringsAsFactors = F)
map.dem = match(dat$patientNum,dat.dem$PatientNum)

# demo
dat$birth_date = dat.dem$Date_of_Birth[map.dem]
dat$death_date = dat.dem$Date_of_Death[map.dem]
dat$gender = dat.dem$Gender[map.dem]
dat$race = dat.dem$Race[map.dem]
dat$marital = dat.dem$Marital_status[map.dem]

write.csv(dat, file = "result/CRC_MAP_dem.csv", 
          row.names = F)
