rm(list=ls())

library(chron)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(survival)
library(foreign)
library(compareGroups)

### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

ic_guapa2<-function(x,y,z)
{
  ic<-paste(x," (",y," to ",z,")",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
  return(pval)
}

pval_guapa2<-function(x)
{
  pval<-ifelse(x<0.00001," < 0.00001",
               ifelse(x<0.001," < 0.001",
                      ifelse(abs(x)<0.01,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                             ifelse(abs(x)<0.1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                                    ifelse(abs(x)<1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),guapa(x))))))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

beta_se_ic_guapa2 <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa2(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa3 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),3)
  ic95a<-round(exp(x-(z*y)),3)
  ic95b<-round(exp(x+(z*y)),3)
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

options(scipen=999)

setwd("C:/Users/Alvaro/Documents/Documentos/Artículos/Montse")
dir.create("./Data")
dir.create("./Outputs")
dir.create("./Outputs/Descriptive")
dir.create("./Outputs/Inflammation")


####################
### DESCRIPTIVES ###
####################

setwd("C:/Users/Alvaro/Documents/Documentos/Artículos/Montse")
load("./Data/predimed_dem_fis01_cc.RData")

xxx<-dat[,c("id","enfneuro","edad","sexo","grup_int","escolar",
            "diabetes","hipercol","hipertg","hta","obesidad","tabaco",
            "p14","actfis","alcohol","apoe4")]
xxx$tot<-1

tab<-NULL
tab<-createTable(compareGroups(enfneuro~.
                               -id-tot,
                               xxx, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                             "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                             "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                 show.n=TRUE, hide.no=0)
tab2<-createTable(compareGroups(tot~.
                                -id-enfneuro,
                                xxx, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                              "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                              "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                  show.n=TRUE, hide.no=0)

tab<-cbind(tab2$descr[,1],tab$descr)
colnames(tab)<-c("All","Non-cases","Cases","P-value","N")
write.table(tab,file="./Outputs/Descriptive/fis01_case_cohort_descriptive.csv",sep=";",col.names=NA)


##############################
### ANALYSES: INFLAMMATION ###
##############################

setwd("C:/Users/Alvaro/Documents/Documentos/Artículos/Montse/Outputs/Inflammation")
dir.create("./splines")
load("C:/Users/Alvaro/Documents/Documentos/Artículos/Montse/predimed_dem_fis01_cc.RData")

# IMPUTATION OF escolar NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))

dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice<-dat_mice[,c("id","escolar")]
dat$escolar<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


# ANALYSES #

vars01<-c("z_log_ldh01","z_c301","z_antitrip01","z_a2m01")
vars02<-c("t_ldh01","t_c301","t_antitrip01","t_a2m01")
ntotal<-7447
z<-qnorm(1-0.05/2)
vars<-c("id","toenfneuro","enfneuro","subcohort","edad","sexo","grup_int","escolar",
        "diabetes","hipercol","hta","tabaco","imc","apoe4")

tab<-NULL
tab2<-NULL

for(i in 1:length(vars01))
  
{  
  vars_def<-append(vars01[i],vars)
  datx<-na.omit(dat[,vars_def])
  nval<-paste(table(datx$enfneuro)[2],"/",table(datx$enfneuro)[1]+table(datx$enfneuro)[2],
              " (",guapa(table(datx$enfneuro)[2]*100/(table(datx$enfneuro)[1]+table(datx$enfneuro)[2])),"%)",sep="")
  
  mod01<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]],
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta01<-exp(summary(mod01)$coefficients[1,1])
  ic95a01<-exp(summary(mod01)$coefficients[1,1]-(z*summary(mod01)$coefficients[1,2]))
  ic95b01<-exp(summary(mod01)$coefficients[1,1]+(z*summary(mod01)$coefficients[1,2]))
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(summary(mod01)$coefficients[1,4])
  
  mod02<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]]
             +edad+as.factor(sexo)+as.factor(grup_int)+as.factor(escolar)
             +as.factor(diabetes)+as.factor(hipercol)
             +as.factor(hta)+as.factor(tabaco)+imc+as.factor(apoe4),
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta02<-exp(summary(mod02)$coefficients[1,1])
  ic95a02<-exp(summary(mod02)$coefficients[1,1]-(z*summary(mod02)$coefficients[1,2]))
  ic95b02<-exp(summary(mod02)$coefficients[1,1]+(z*summary(mod02)$coefficients[1,2]))
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(summary(mod02)$coefficients[1,4])
  
  datx[,vars02[i]]<-cut(datx[,vars01[i]],quantile(datx[,vars01[i]],probs=seq(0,1,1/3),na.rm=TRUE),labels=FALSE)
  datx[,vars02[i]]<-with(datx,ifelse(datx[,vars02[i]]==1,0,
                                     ifelse(datx[,vars02[i]]==2,NA,
                                            ifelse(datx[,vars02[i]]==3,1,NA))))
  datx<-subset2(datx,"!is.na(datx[,vars02[i]])")
  
  mod04<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars02[i]]
             +edad+as.factor(sexo)+as.factor(grup_int)+as.factor(escolar)
             +as.factor(diabetes)+as.factor(hipercol)
             +as.factor(hta)+as.factor(tabaco)+imc+as.factor(apoe4),
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta04<-exp(summary(mod04)$coefficients[1,1])
  ic95a04<-exp(summary(mod04)$coefficients[1,1]-(z*summary(mod04)$coefficients[1,2]))
  ic95b04<-exp(summary(mod04)$coefficients[1,1]+(z*summary(mod04)$coefficients[1,2]))
  coef04<-ic_guapa2(guapa(beta04),guapa(ic95a04),guapa(ic95b04))
  pval04<-pval_guapa(summary(mod04)$coefficients[1,4])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef04,pval04))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("HR (raw)","pval (raw)","HR (mod2)","pval (mod2)","HR (tert-mod2)","pval (tert-mod2)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./case_cohort_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./case_cohort_forestplots.csv",sep=";",col.names=NA)


dat<-read.csv2("./case_cohort_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,
             2,2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-with(dat,ifelse(X=="z_log_ldh01","01.\nLactate\ndehydrogenase",
                         ifelse(X=="z_c301","02.\nComplement\ncomponent 3",
                                ifelse(X=="z_antitrip01","03.\nAlpha-1\nantitripsin",
                                       ifelse(X=="z_a2m01","04.\nAlpha-2\nmacroglobulin",NA)))))
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./case_cohort_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(2.5)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./case_cohort_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + ylab("Risk of dementia\nper 1 SD higher levels of biomarkers at baseline\n(hazard ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.5,0.75,1,1.5,2)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=11, scales="free_y") +
    coord_flip() +
    theme_minimal() +
    theme(legend.position="bottom",
          #panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=11),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11),
          strip.text=element_text(size=6.4))
  
  namefile<-paste("./case_cohort_",vars01[i],".tiff",sep="")
  ggsave(filename=namefile, units="px", width=9000, height=5000, dpi=1200, bg="transparent")
  figure
}


### CAMBIOS A 1 ANY ###

load("C:/Users/Alvaro/Documents/Documentos/Artículos/Montse/Data/predimed_dem_fis01_cc.RData")

# IMPUTATION OF LIPID PROFILE NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))

dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice<-dat_mice[,c("id","escolar")]
dat$escolar<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


# ANALYSES #

vars01<-c("z_log_ldh01","z_c301","z_antitrip01","z_a2m01")
vars02<-c("z_log_ldh03","z_c303","z_antitrip03","z_a2m03")

ntotal<-7447
z<-qnorm(1-0.05/2)
vars<-c("id","toenfneuro","enfneuro","subcohort","edad","sexo","grup_int","escolar",
        "diabetes","hipercol","hta","tabaco","imc","apoe4")

tab<-NULL
tab2<-NULL

for(i in 1:length(vars01))
  
{  
  vars_def<-append(vars01[i],vars)
  vars_def<-append(vars02[i],vars_def)
  datx<-na.omit(dat[,vars_def])
  nval<-paste(table(datx$enfneuro)[2],"/",table(datx$enfneuro)[1]+table(datx$enfneuro)[2],
              " (",guapa(table(datx$enfneuro)[2]*100/(table(datx$enfneuro)[1]+table(datx$enfneuro)[2])),"%)",sep="")
  
  mod01<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]]+datx[,vars02[i]],
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta01<-exp(summary(mod01)$coefficients[1,1])
  ic95a01<-exp(summary(mod01)$coefficients[1,1]-(z*summary(mod01)$coefficients[1,2]))
  ic95b01<-exp(summary(mod01)$coefficients[1,1]+(z*summary(mod01)$coefficients[1,2]))
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(summary(mod01)$coefficients[1,4])
  
  mod02<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]]+datx[,vars02[i]]
             +edad+as.factor(sexo)+as.factor(grup_int)+as.factor(escolar)
             +as.factor(diabetes)+as.factor(hipercol)
             +as.factor(hta)+as.factor(tabaco)+imc+as.factor(apoe4),
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta02<-exp(summary(mod02)$coefficients[1,1])
  ic95a02<-exp(summary(mod02)$coefficients[1,1]-(z*summary(mod02)$coefficients[1,2]))
  ic95b02<-exp(summary(mod02)$coefficients[1,1]+(z*summary(mod02)$coefficients[1,2]))
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(summary(mod02)$coefficients[1,4])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("HR (raw)","pval (raw)","HR (mod2)","pval (mod2)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./case_cohort_1year_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./case_cohort_1year_forestplots.csv",sep=";",col.names=NA)


dat<-read.csv2("./case_cohort_1year_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,
             2,2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-with(dat,ifelse(X=="z_log_ldh01","01.\nLactate\ndehydrogenase",
                         ifelse(X=="z_c301","02.\nComplement\ncomponent 3",
                                ifelse(X=="z_antitrip01","03.\nAlpha-1\nantitripsin",
                                       ifelse(X=="z_a2m01","04.\nAlpha-2\nmacroglobulin",NA)))))
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./case_cohort_1year_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(3)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./case_cohort_1year_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + ylab("Risk of dementia\nper 1 SD increases in 1-year differences in biomarkers\n(hazard ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.25,0.5,1,2,4)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=11, scales="free_y") +
    coord_flip() +
    theme_minimal() +
    theme(legend.position="bottom",
          #panel.spacing = unit(1, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_blank(),
          axis.text.x=element_text(size=11),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11),
          strip.text=element_text(size=6.4))
  
  namefile<-paste("./case_cohort_1year_",vars01[i],".tiff",sep="")
  ggsave(filename=namefile, units="px", width=9000, height=5000, dpi=1200, bg="transparent")
  figure
}

