
dat.trt = read.csv( "data/paper surg/CRC_surg_stage_bmi.csv",
                    stringsAsFactors = F)
dat.bmi.nlp =read.csv("data/paper surg/RCTsurg_BMI_clean_NLP.csv",
                               stringsAsFactors = F)
dat.bmi.str = read.csv("data/paper surg/RCTsurg_BMI_clean_Structure.csv",
                               stringsAsFactors = F)

dat.trt$bmi.str = NA
dat.trt$bmi.nlp = NA

map.pos = match(dat.bmi.str$PatientNum, dat.trt$patientNum)
best.pos = !duplicated(dat.bmi.str$PatientNum)
dat.trt$bmi.str[map.pos[best.pos]] = dat.bmi.str$bmi.med365[best.pos]


map.pos = match(dat.bmi.nlp$PatientNum, dat.trt$patientNum)
best.pos = !duplicated(dat.bmi.nlp$PatientNum)
dat.trt$bmi.nlp[map.pos[best.pos]] = dat.bmi.nlp$bmi.med365[best.pos]


dat.bmi.both = dat.trt[!is.na(dat.trt$bmi.str) & !is.na(dat.trt$bmi.nlp),]
dat.bmi.both.emu = dat.bmi.both[dat.bmi.both$patientNum %in% dat$patientNum,]

plot(dat.bmi.both$bmi.nlp, dat.bmi.both$bmi.str)
plot(dat.bmi.both.emu$bmi.nlp, dat.bmi.both.emu$bmi.str)

library(ggplot2)

ggdf = data.frame(Structured = dat.bmi.both$bmi.str, 
                  NLP = dat.bmi.both$bmi.nlp)
p1 = ggplot(ggdf, aes(x=Structured, y = NLP))
p1 = p1 + geom_point()
p1 = p1 + geom_abline(slope = 1, intercept = 0)
p1 = p1 + labs(title = "All Patients")
p1 = p1 + geom_abline(intercept = 30, slope = 0, linetype = "dashed")
p1 = p1 + geom_vline(xintercept = 30, linetype = "dashed")
p1 = p1 + geom_abline(intercept = 25, slope = 0, linetype = "dashed")
p1 = p1 + geom_vline(xintercept = 25, linetype = "dashed")


ggdf = data.frame(Structured = dat.bmi.both.emu$bmi.str, 
                  NLP = dat.bmi.both.emu$bmi.nlp)
p2 = ggplot(ggdf, aes(x=Structured, y = NLP))
p2 = p2 + geom_point()
p2 = p2 + geom_abline(slope = 1, intercept = 0)
p2 = p2 + labs(title = "Emulation Cohort Patients")
p2 = p2 + geom_abline(intercept = 30, slope = 0, linetype = "dashed")
p2 = p2 + geom_vline(xintercept = 30, linetype = "dashed")
p2 = p2 + geom_abline(intercept = 25, slope = 0, linetype = "dashed")
p2 = p2 + geom_vline(xintercept = 25, linetype = "dashed")
# Missclassification: 
table(dat.bmi.both.emu$bmi.nlp>30, dat.bmi.both.emu$bmi.str>30)

library(gridExtra)
library(grid)

png("figures/paper surg/BMI_concord.png", 
    width =2400, height = 1200, res = 300)
grid.arrange(arrangeGrob(p1,p2,
                         layout_matrix = matrix(1:2,1,2),
                         widths = c(1,1)
),
ncol = 1
)
dev.off()