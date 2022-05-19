# local.path = "D:/Local folders for dropbox/"
local.path = "E:/Local research data/"

dat =read.csv("data/paper surg/CRC_surg_features.csv",
              stringsAsFactors = F)
dat.CPT = read.csv("data/paper surg/CRC_surg_CPT.csv", 
                   stringsAsFactors = F)
LAC.gen = "C44204"
OC.gen = "C44140"
LAC.right = "C44205"
OC.right = "C44160"
LAC.sigmoid = c("C44207", "C44208")
OC.sigmoid = c("C44145", "C44146")
CPT.comon = c(LAC.gen, OC.gen, LAC.right, OC.right, 
              LAC.sigmoid, OC.sigmoid)
dat$OC.CPT = dat.CPT$open.filtered.code[match(dat$patientNum, dat.CPT$patientNum)]
dat$LAC.CPT = dat.CPT$lap.filtered.code[match(dat$patientNum, dat.CPT$patientNum)]
dat$CPT = dat$OC.CPT
dat$CPT[dat$lap.CPT.filtered] = dat$LAC.CPT[dat$lap.CPT.filtered]


date.var = c("birth_date","death_date","last.date","trt_date", "first_CRC_dx")
for (tmpvar in date.var)
{
  dat[,tmpvar] =as.Date(dat[,tmpvar])
}
dat = dat[(dat$trt_date <= as.Date("2017-12-31"))&
            (dat$trt_date >= as.Date("2006-01-01")),]
dat = dat[as.numeric(dat$last.date - dat$trt_date, units = "days")>30, 
]

gap = as.numeric(dat$trt_date - dat$first_CRC_dx, units = "days")
table( gap>=0 & gap <= 90)
dat = dat[gap>=0 & gap <= 90, ]

dat.rct = read.csv(paste(local.path,"Merck/RCT data/NCT00002575-Colectomy/master.csv", sep=''),
                   stringsAsFactors = F)
dat.rct = dat.rct[!is.na(dat.rct$fu_yrs),]
dat.rct = dat.rct[dat.rct$stage!=4,]
dat.rct = dat.rct[!is.na(dat.rct$bow_adh),]

# Apply some eligibility criteria
dat = dat[as.numeric(format(dat$trt_date, "%Y")) %in% 2006:2017, ]

table(dat$crohn_new,dat$lap.CPT.filtered) # 43
dat = dat[!dat$crohn_new,] # Crohn's disease rare 14,10
table(dat$cuc_new,dat$lap.CPT.filtered) # 40
dat = dat[!dat$cuc_new,] # rare 5,1
table(dat$CRC_seg_surgday,dat$lap.CPT.filtered)
dat = dat[!(dat$CRC_seg_surgday %in% c("multiple", "rectum","unknown","transverse colon")), ] # rare in LAC 4,9
table(dat$CRC_seg_surgday, dat$CPT)
table(!(dat$CPT %in% CPT.comon),dat$lap.CPT.filtered)
dat = dat[dat$CPT %in% CPT.comon ,]
table((dat$CRC_seg_surgday=="Right" & 
         dat$CPT %in% c(OC.sigmoid, LAC.sigmoid))
      | (dat$CRC_seg_surgday=="sigmoid" & 
           dat$CPT %in% c(OC.right, LAC.right))
      | (dat$CRC_seg_surgday=="Left" & 
           dat$CPT %in% c(OC.sigmoid, LAC.sigmoid,
                          OC.right, LAC.right)),
      dat$lap.CPT.filtered)
dat = dat[!(dat$CRC_seg_surgday=="Right" & 
              dat$CPT %in% c(OC.sigmoid, LAC.sigmoid)), ]
dat = dat[!(dat$CRC_seg_surgday=="sigmoid" & 
              dat$CPT %in% c(OC.right, LAC.right)), ]
dat = dat[!(dat$CRC_seg_surgday=="Left" & 
              dat$CPT %in% c(OC.sigmoid, LAC.sigmoid,
                             OC.right, LAC.right)), ]
table(dat$obstruct_new, dat$lap.CPT.filtered) # 537
dat = dat[!dat$obstruct_new,] # rare in LAC 13
table(dat$stage.wICD, dat$lap.CPT.filtered)
dat = dat[-which(dat$stage.wICD==4),]
sum(is.na(dat$stage.wICD))
dat = dat[!is.na(dat$stage.wICD),]
table(dat$prior_cancer_new, dat$lap.CPT.filtered)
dat = dat[!dat$prior_cancer_new,]

table(dat$adhesion, dat$lap.CPT.filtered)
adhesion = as.numeric(dat$adhesion_new)
seg = factor(dat$CRC_seg_surgday, levels =  c("Right", "Left", "sigmoid"))
right.gen = as.numeric(dat$CRC_seg_surgday=="Right" & 
                         dat$CPT %in% c(OC.gen, LAC.gen))
sigmoid.gen = as.numeric(dat$CRC_seg_surgday=="sigmoid" & 
                           dat$CPT %in% c(OC.gen, LAC.gen))
trt = as.numeric(dat$lap.CPT.filtered)
stage = factor(dat$stage.wICD, levels = c(0:3),
               exclude = NULL,
               labels = c("0-1","0-1",2:3)) # Merge O-I, O rare, 2,0

# Follow-up duration for MAP
fu = as.numeric(as.Date(dat$last.date)-as.Date(dat$first.date), units = "days")/365.2425

# Pre treatment follow-up
pre.fu = as.numeric(as.Date(dat$last.date)-dat$trt_date, units = "days")/365.2425

library(ggplot2)

ggdf = data.frame(x = fu[trt==1])
p1 <- ggplot(ggdf, aes(x=x))+
  geom_histogram(alpha=0.5, position="identity", color="black", bins = 20)#+
p1 <- p1 + labs(title = "LAC arm")
p1 <- p1 + scale_x_continuous("Total EHR Follow-up (Years)",seq(0,40,5))
p1 <- p1 + ylim(0,70)

ggdf = data.frame(x = fu[trt==0])
p2 <- ggplot(ggdf, aes(x=x))+
  geom_histogram(alpha=0.5, position="identity", color="black", bins = 20)#+
p2 <- p2 + labs(title = "OC arm")
p2 <- p2 + scale_x_continuous("Total EHR Follow-up (Years)",seq(0,40,5))
p2 <- p2 + ylim(0,70)


ggdf = data.frame(x = pre.fu[trt==1])
p3 <- ggplot(ggdf, aes(x=x))+
  geom_histogram(alpha=0.5, position="identity", color="black", bins = 10)#+
p3 <- p3 + labs(title = "LAC arm")
p3 <- p3 + scale_x_continuous("Pre-colectomy Follow-up (Years)",seq(0,40,5))
p3 <- p3 + ylim(0,130)

ggdf = data.frame(x = pre.fu[trt==0])
p4 <- ggplot(ggdf, aes(x=x))+
  geom_histogram(alpha=0.5, position="identity", color="black", bins = 10)#+
p4 <- p4 + labs(title = "OC arm")
p4 <- p4 + scale_x_continuous("Pre-colectomy Follow-up (Years)",seq(0,40,5))
p4 <- p4 + ylim(0,130)

library(gridExtra)
library(grid)

png("figures/paper surg/follow_up_years.png", 
    width = 1800, height = 1800, res = 300)
grid.arrange(arrangeGrob(p1,p2,
                         p3, p4,
                         layout_matrix = t(matrix(1:4,2,2)),
                         widths = c(1,1)
),
ncol = 1
)
dev.off()