setwd("E:/Dropbox/HDLBP (1)/")

dat<-read.delim("Xenograft/xenograft_all.txt", stringsAsFactors = F, dec = ",")

library(ggplot2)

dat$status<-factor(dat$status, levels=c("WT", "KO"))
dat$time<-factor(dat$time, levels=c("7", "14", "21"))


ggplot(dat, aes(time, log2(Volume), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor volume")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(time, (Volume), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("tumor volume")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(time, log2(Width), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor width")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))
ggplot(dat, aes(time, log2(Height), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor height")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))
ggplot(dat, aes(time, (Width), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("tumor width")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(time, log2(Length), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor length")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))
ggplot(dat, aes(time, (Length), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.5))+xlab("time (days)")+ylab("log2 tumor length")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(status, (tumor_weight_mg), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.8))+xlab("")+ylab("tumor weight (mg)")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))

ggplot(dat, aes(status, log2(tumor_weight_mg), colour=status))+geom_boxplot()+
  geom_point( position = position_jitterdodge(jitter.width = 0.8))+xlab("")+ylab("log2 tumor weight (mg)")+
  scale_colour_manual(values=c("dodgerblue3", "orange3"))


t.test(dat[dat$status=="WT" & dat$time=="7", "Volume"], dat[dat$status=="KO" & dat$time=="7", "Volume"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Volume"], dat[dat$status=="KO" & dat$time=="14", "Volume"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Volume"], dat[dat$status=="KO" & dat$time=="21", "Volume"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Length"], dat[dat$status=="KO" & dat$time=="7", "Length"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Length"], dat[dat$status=="KO" & dat$time=="14", "Length"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Length"], dat[dat$status=="KO" & dat$time=="21", "Length"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Width"], dat[dat$status=="KO" & dat$time=="7", "Width"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Width"], dat[dat$status=="KO" & dat$time=="14", "Width"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Width"], dat[dat$status=="KO" & dat$time=="21", "Width"])

t.test(dat[dat$status=="WT" & dat$time=="7", "Height"], dat[dat$status=="KO" & dat$time=="7", "Height"])
t.test(dat[dat$status=="WT" & dat$time=="14", "Height"], dat[dat$status=="KO" & dat$time=="14", "Height"])
t.test(dat[dat$status=="WT" & dat$time=="21", "Height"], dat[dat$status=="KO" & dat$time=="21", "Height"])


t.test(dat[dat$status=="WT" & dat$time=="21", "tumor_weight_mg"], dat[dat$status=="KO" & dat$time=="21", "tumor_weight_mg"])


length(na.omit(dat[dat$status=="WT" & dat$time=="7", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="7", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="14", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="14", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="21", "Volume"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="21", "Volume"]))
length(na.omit(dat[dat$status=="WT" & dat$time=="21", "tumor_weight_mg"]))
length(na.omit(dat[dat$status=="KO" & dat$time=="21", "tumor_weight_mg"]))




anova_two_way<-aov(Volume~time+status, data=dat)
summary(anova_two_way)

anova_two_way<-aov(Length~time+status, data=dat)
summary(anova_two_way)

anova_two_way<-aov(Width~time+status, data=dat)
summary(anova_two_way)

anova_two_way<-aov(Height~time+status, data=dat)
summary(anova_two_way)

anova_one_way<-aov(tumor_weight_mg~status, data=dat)
summary(anova_one_way)


anova_two_way_int<-aov(Volume~time*status, data=dat)
summary(anova_two_way_int)

TukeyHSD(anova_two_way_int, which = c("time","status"))

dat$status_time<-paste0(dat$status,"_",dat$time)

anova_two_way_int<-aov(Volume~time+status_time, data=dat)
summary(anova_two_way_int)

TukeyHSD(anova_two_way_int, which = "status_time")


