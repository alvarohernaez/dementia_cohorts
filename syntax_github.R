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


header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

options(scipen=999)


setwd(".../PREDIMED_dementia")
dir.create("./Data")
dir.create("./Outputs")
dir.create("./Outputs/Descriptive")
dir.create("./Outputs/Thrombosis")
dir.create("./Outputs/Adiposity")
dir.create("./Outputs/Lipids")
dir.create("./Outputs/Inflammation")


####################
### DESCRIPTIVES ###
####################

### FIS01 - CASE COHORT ###
###########################

### DESCRIPTIVE ###

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


### SELECTION BIAS ###

sel<-as.data.frame(dat[,c("id")])
sel$sel<-1
names(sel)<-c("id","sel")

load("./Data/predimed_dem_total.RData")
dat<-merge2(dat,sel,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$sel<-with(dat,ifelse(is.na(sel),0,sel))

dat<-dat[,c("id","edad","sexo","grup_int","escolar",
            "diabetes","hipercol","hipertg","hta","obesidad","tabaco",
            "p14","actfis","alcohol","apoe4","sel")]

tab<-NULL
tab<-createTable(compareGroups(sel~.
                               -id,
                               dat, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                             "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                             "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                 show.n=TRUE, hide.no=0)

tab<-cbind(tab$descr)
colnames(tab)<-c("Non-selected","Selected","P-value","N")
write.table(tab,file="./Outputs/Descriptive/fis01_case_cohort_selection_bias.csv",sep=";",col.names=NA)


### FIS02 - NESTED CASE CONTROL ###
###################################

### DESCRIPTIVE ###

load("./Data/predimed_dem_fis02_ncc.RData")

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
write.table(tab,file="./Outputs/Descriptive/fis02_nested_case_control_descriptive.csv",sep=";",col.names=NA)


### SELECTION BIAS ###

sel<-as.data.frame(dat[,c("id")])
sel$sel<-1
names(sel)<-c("id","sel")

load("./Data/predimed_dem_total.RData")
dat<-merge2(dat,sel,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$sel<-with(dat,ifelse(is.na(sel),0,sel))

dat<-dat[,c("id","edad","sexo","grup_int","escolar",
            "diabetes","hipercol","hipertg","hta","obesidad","tabaco",
            "p14","actfis","alcohol","apoe4","sel")]

tab<-NULL
tab<-createTable(compareGroups(sel~.
                               -id,
                               dat, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                             "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                             "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                 show.n=TRUE, hide.no=0)

tab<-cbind(tab$descr)
colnames(tab)<-c("Non-selected","Selected","P-value","N")
write.table(tab,file="./Outputs/Descriptive/fis02_nested_case_control_selection_bias.csv",sep=";",col.names=NA)


### PLATELETS ###
#################

### DESCRIPTIVE ###

load("./Data/predimed_dem_total.RData")
#plot(compareGroups(~logtg,data=dat))

dat<-dat[!is.na(dat$plaquetas),] # n = 4393

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
write.table(tab,file="./Outputs/Descriptive/plaquetas_descriptive.csv",sep=";",col.names=NA)


### SELECTION BIAS ###

sel<-as.data.frame(dat[,c("id")])
sel$sel<-1
names(sel)<-c("id","sel")

load("./Data/predimed_dem_total.RData")
dat<-merge2(dat,sel,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$sel<-with(dat,ifelse(is.na(sel),0,sel))

dat<-dat[,c("id","edad","sexo","grup_int","escolar",
            "diabetes","hipercol","hipertg","hta","obesidad","tabaco",
            "p14","actfis","alcohol","apoe4","sel")]

tab<-NULL
tab<-createTable(compareGroups(sel~.
                               -id,
                               dat, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                             "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                             "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                 show.n=TRUE, hide.no=0)

tab<-cbind(tab$descr)
colnames(tab)<-c("Non-selected","Selected","P-value","N")
write.table(tab,file="./Outputs/Descriptive/plaquetas_selection_bias.csv",sep=";",col.names=NA)



### LIPID PROFILE ###
#####################

### DESCRIPTIVE ###

load("./Data/predimed_dem_total.RData")

dat<-dat[!is.na(dat$tc) | !is.na(dat$hdlc) | !is.na(dat$logtg),] # n = 7079

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
write.table(tab,file="./Outputs/Descriptive/lipids_descriptive.csv",sep=";",col.names=NA)


### SELECTION BIAS ###

sel<-as.data.frame(dat[,c("id")])
sel$sel<-1
names(sel)<-c("id","sel")

load("./Data/predimed_dem_total.RData")
dat<-merge2(dat,sel,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$sel<-with(dat,ifelse(is.na(sel),0,sel))

dat<-dat[,c("id","edad","sexo","grup_int","escolar",
            "diabetes","hipercol","hipertg","hta","obesidad","tabaco",
            "p14","actfis","alcohol","apoe4","sel")]

tab<-NULL
tab<-createTable(compareGroups(sel~.
                               -id,
                               dat, method=c("sexo"=3,"grup_int"=3,"apoe4"=3,"escolar"=3,
                                             "diabetes"=3,"hipercol"=3,"hipertg"=3,"hta"=3,
                                             "tabaco"=3,"obesidad"=3,"actfis"=2,"alcohol"=2)),
                 show.n=TRUE, hide.no=0)

tab<-cbind(tab$descr)
colnames(tab)<-c("Non-selected","Selected","P-value","N")
write.table(tab,file="./Outputs/Descriptive/lipids_selection_bias.csv",sep=";",col.names=NA)



### ALL PARTICIPANTS: IMC, CINTURA, MEDICATION ###
##################################################

### DESCRIPTIVE ###

load("./Data/predimed_dem_total.RData")

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
write.table(tab,file="./Outputs/Descriptive/all_participants_descriptive.csv",sep=";",col.names=NA)



############################
### ANALYSES: THROMBOSIS ###
############################

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Outputs/Thrombosis")
dir.create("./splines")


### NESTED CASE CONTROL: CONDITIONAL LOGISTIC REGRESSIONS ###
#############################################################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis02_ncc.RData")

# PREPARATION OF COVARIATES / IMPUTATION OF ESCOLAR AND PLATELET NAs #

round(sum(is.na(dat$plaquetas))/dim(dat)[1]*100,0) #39% NAs in plaquetas
sum(dat$escolar==9)

dat_mice<-dat[,c("id","enfneuro",
                 "sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc","plaquetas","antiplaq","anticoagoral")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
round(sum(is.na(dat_mice$escolar))/dim(dat_mice)[1]*100,1) #0.4% NAs in escolar

dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice<-dat_mice[,c("id","escolar","plaquetas")]
dat$escolar<-NULL
dat$plaquetas<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


# ANALYSES #

vars01<-c("z_log_pselectin","z_log_fibrinogen","z_pai1","z_factorv","z_factorvii","z_log_factorviii")
vars02<-c("log_pselectin","log_fibrinogen","pai1","factorv","factorvii","log_factorviii")
vars03<-c("P-selectin, ng/mL (log-transformed)","Fibrinogen, mg/dL (log-transformed)",
          "PAI-1, ng/mL","Factor V, mg/dL","Factor VII, mg/dL","Factor VIII, mg/dL (log-transformed)")

z<-qnorm(1-0.05/2)
closest<-function(xv,sv)
{
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]
}

tab<-NULL
tab2<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars01[i]])")
  
  aaa<-datx[,vars01[i]]
  
  mod01<-clogit(enfneuro~datx[,vars01[i]]+strata(match),
                data=datx)
  beta01<-exp(summary(mod01)$coefficients[1,1])
  ic95a01<-exp(summary(mod01)$coefficients[1,1]-(z*summary(mod01)$coefficients[1,3]))
  ic95b01<-exp(summary(mod01)$coefficients[1,1]+(z*summary(mod01)$coefficients[1,3]))
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(summary(mod01)$coefficients[1,5])
  
  mod02<-clogit(enfneuro~datx[,vars01[i]]+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                +as.factor(apoe4),
                data=datx)
  beta02<-exp(summary(mod02)$coefficients[1,1])
  ic95a02<-exp(summary(mod02)$coefficients[1,1]-(z*summary(mod02)$coefficients[1,3]))
  ic95b02<-exp(summary(mod02)$coefficients[1,1]+(z*summary(mod02)$coefficients[1,3]))
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(summary(mod02)$coefficients[1,5])
  
  mod03<-clogit(enfneuro~datx[,vars01[i]]+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                +plaquetas+as.factor(antiplaq)+as.factor(anticoagoral)+as.factor(apoe4),
                data=datx)
  beta03<-exp(summary(mod03)$coefficients[1,1])
  ic95a03<-exp(summary(mod03)$coefficients[1,1]-(z*summary(mod03)$coefficients[1,3]))
  ic95b03<-exp(summary(mod03)$coefficients[1,1]+(z*summary(mod03)$coefficients[1,3]))
  coef03<-ic_guapa2(guapa(beta03),guapa(ic95a03),guapa(ic95b03))
  pval03<-pval_guapa(summary(mod03)$coefficients[1,5])
  
  mod03_nl<-clogit(enfneuro~bs(aaa,degree=3)+strata(match)
                   +as.factor(grup_int)+as.factor(escolar)
                   +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                   +plaquetas+as.factor(antiplaq)+as.factor(anticoagoral)+as.factor(apoe4),
                   data=datx)
  p_nonlin<-pval_guapa(lrtest(mod03,mod03_nl)[2,5])
  
  aaa<-datx[,vars02[i]]
  mod03<-clogit(enfneuro~aaa+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                +plaquetas+as.factor(antiplaq)+as.factor(anticoagoral)+as.factor(apoe4),
                data=datx)
  mod03_nl<-clogit(enfneuro~bs(aaa,degree=3)+strata(match)
                   +as.factor(grup_int)+as.factor(escolar)
                   +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                   +plaquetas+as.factor(antiplaq)+as.factor(anticoagoral)+as.factor(apoe4),
                   data=datx)
  
  ptemp<-termplot(mod03_nl,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,median(aaa,na.rm=TRUE))[1]
  center<-with(temp, y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./splines/",vars02[i],".jpg",sep="")
  labely<-c("Neurodegenerative disease risk (odds ratio)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  plot.data<-subset2(plot.data,"plot.data$lci>=0.1 & plot.data$uci<=10")
  
  p_lin2<-pval_guapa(summary(mod03)$coefficients[1,5])
  p_nonlin2<-pval_guapa(lrtest(mod03,mod03_nl)[2,5])
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,
             "\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(data=plot.data, aes_string(x=plot.data$x, y=plot.data$y)) + 
    geom_ribbon(aes_string(ymin=plot.data$lci, ymax=plot.data$uci), alpha=0.25, fill="Black") +
    geom_line(aes_string(x=plot.data$x, y=plot.data$y), color='Black') +
    geom_hline(yintercept=1, linetype=2) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    labs(x=vars03[i],y=labely) +
    annotate("text", x=max(plot.data$x,na.rm=TRUE)*0.98, y=max(plot.data$uci,na.rm=TRUE), label=leg, vjust=1, hjust=1, size=6) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  
  ggsave(filename=name, dpi=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,
                       coef03,pval03,p_nonlin))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02,
                         beta03,ic95a03,ic95b03))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (mod2)","pval (mod2)",
                 "OR (mod3)","pval (mod3)","pval-nonlin (adj)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./nested_case_control_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./nested_case_control_forestplots.csv",sep=";",col.names=NA)


### GENERATION OF FOREST PLOT ###

dat<-read.csv2("./nested_case_control_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)]),
                         as.matrix(dat[c(1,8,9,10)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,1,1,
             2,2,2,2,2,2,
             3,3,3,3,3,3)
dat$model<-with(dat,ifelse(model==1,"03. Not adjusted",
                           ifelse(model==2,"02. Adjusted",
                                  ifelse(model==3,"01. +plaq+med",NA))))
dat$bmk<-c("01.\nP-selectin","02.\nFibrinogen","03.\nPAI-1",
           "04.\nFactor V","05.\nFactor VII","06.\nFactor VIII",
           "01.\nP-selectin","02.\nFibrinogen","03.\nPAI-1",
           "04.\nFactor V","05.\nFactor VII","06.\nFactor VIII",
           "01.\nP-selectin","02.\nFibrinogen","03.\nPAI-1",
           "04.\nFactor V","05.\nFactor VII","06.\nFactor VIII")
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./nested_case_control_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(3.5)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./nested_case_control_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5.5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + 
    ylab("Risk of dementia\nper 1 SD higher levels of biomarkers at baseline\n(odds ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.25,0.5,1,2,4)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=13, scales="free_y") +
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
          axis.title.y=element_text(size=20),
          strip.text=element_text(size=10))
  
  ggsave("./nested_case_control_forestplot.tiff", units="px", width=9000, height=8000, dpi=1200, bg="transparent")
  figure
}



### SURVIVAL ANALYSES ###
#########################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_total.RData")

# IMPUTATION OF escolar AND apoe4 NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-with(dat_mice,ifelse(apoe4==9,NA,apoe4))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
round(sum(is.na(dat_mice$escolar))/dim(dat_mice)[1]*100,1) #1.8% NAs in escolar
round(sum(is.na(dat_mice$apoe4))/dim(dat_mice)[1]*100,1) #3% NAs in escolar

dat_mice[] <- lapply(dat_mice, function(x) if (is.factor(x) & !identical(x, dat_mice$escolar) & !identical(x, dat_mice$apoe4)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
dat_mice<-dat_mice[,c("id","escolar","apoe4")]
dat$escolar<-NULL
dat$apoe4<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)


dat$adjx<-1

z<-qnorm(1-0.05/2)

vars01<-c("z_plaquetas","antiplaq","anticoagoral")
vars02<-c("antiplaq","anticoagoral","adjx")
vars03<-c("anticoagoral","adjx","antiplaq")

tab<-NULL
tab2<-NULL


for(i in 1:length(vars01))
  
{ 
  participants<-length(which(!is.na(dat[,vars01[i]])))
  dat2<-dat[!is.na(dat[,vars01[i]]),]
  
  mod01<-coxph(Surv(toenfneuro, enfneuro)~dat2[,vars01[i]]+cluster(idcluster2),
               na.action=na.exclude, method="breslow", data=dat2)
  beta01<-intervals(mod01)[1,1]
  ic95a01<-intervals(mod01)[1,2]
  ic95b01<-intervals(mod01)[1,3]
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(intervals(mod01)[1,4])
  
  mod02<-coxph(Surv(toenfneuro, enfneuro)~dat2[,vars01[i]]+cluster(idcluster2)
               +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
               +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
               +dat2[,vars02[i]]+dat2[,vars03[i]]+as.factor(apoe4),
               na.action=na.exclude, method="breslow", data=dat2)
  beta02<-intervals(mod02)[1,1]
  ic95a02<-intervals(mod02)[1,2]
  ic95b02<-intervals(mod02)[1,3]
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(intervals(mod02)[1,4])
  
  #aaa<-dat2[,vars01[i]]
  #mfit<-coxph(Surv(toenfneuro, enfneuro)~pspline(aaa, df=4)
  #            +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
  #            +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
  #            +dat2[,vars02[i]]+dat2[,vars03[i]]+as.factor(apoe4),
  #            na.action=na.exclude, method="breslow",data=dat2)
  #p_nonlin<-pval_guapa(lrtest(mod02,mfit)[2,5])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (adj)","pval (adj)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./platelets_med_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./platelets_med_forestplots.csv",sep=";",col.names=NA)


### GENERATION OF FOREST PLOT ###

dat<-read.csv2("./platelets_med_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,
             2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-c("01.\nPlatelets\n(+1 SD)","02.\nAntiplatelet\ndrug use\n(vs. non-use)",
           "03.\nVit. K antagonist\ndrug use\n(vs. non.use)",
           "01.\nPlatelets\n(+1 SD)","02.\nAntiplatelet\ndrug use\n(vs. non-use)",
           "03.\nVit. K antagonist\ndrug use\n(vs. non.use)")
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./platelets_med_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(4.5)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./platelets_med_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5.5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + 
    ylab("Risk of dementia\n(hazard ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.25,0.5,1,2,4)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=13, scales="free_y") +
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
          axis.title.y=element_text(size=20),
          strip.text=element_text(size=6.5))
  
  ggsave("./platelets_med_forestplot.tiff", units="px", width=9000, height=4000, dpi=1200, bg="transparent")
  figure
}



###########################
### ANALYSES: ADIPOSITY ###
###########################

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Outputs/Adiposity")
dir.create("./splines")


### NESTED CASE CONTROL: CONDITIONAL LOGISTIC REGRESSIONS ###
#############################################################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis02_ncc.RData")

# PREPARATION OF COVARIATES / IMPUTATION OF ESCOLAR NAs #

sum(dat$escolar==9)

dat_mice<-dat[,c("id","enfneuro",
                 "sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))

dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice<-dat_mice[,c("id","escolar")]
dat$escolar<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


### BIOMARCADORES ###
#####################

vars01<-c("z_leptin","z_ghrelin","z_resistin","z_adipsin","z_log_fgf21")
vars02<-c("leptin","ghrelin","resistin","adipsin","log_fgf21")
vars03<-c("Leptin, ng/mL","Ghrelin, pg/mL","Resistin, ng/mL","Adipsin, ng/mL",
          "Fibroblast growth factor-21, pg/mL (log-transformed)")

z<-qnorm(1-0.05/2)

tab<-NULL
tab2<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"!is.na(dat[,vars01[i]])")
  
  aaa<-datx[,vars01[i]]
  
  mod01<-clogit(enfneuro~datx[,vars01[i]]+strata(match),
                data=datx)
  beta01<-exp(summary(mod01)$coefficients[1,1])
  ic95a01<-exp(summary(mod01)$coefficients[1,1]-(z*summary(mod01)$coefficients[1,3]))
  ic95b01<-exp(summary(mod01)$coefficients[1,1]+(z*summary(mod01)$coefficients[1,3]))
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(summary(mod01)$coefficients[1,5])
  
  mod02<-clogit(enfneuro~datx[,vars01[i]]+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)
                +as.factor(apoe4),
                data=datx)
  beta02<-exp(summary(mod02)$coefficients[1,1])
  ic95a02<-exp(summary(mod02)$coefficients[1,1]-(z*summary(mod02)$coefficients[1,3]))
  ic95b02<-exp(summary(mod02)$coefficients[1,1]+(z*summary(mod02)$coefficients[1,3]))
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(summary(mod02)$coefficients[1,5])
  
  mod03<-clogit(enfneuro~datx[,vars01[i]]+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)
                +imc+as.factor(apoe4),
                data=datx)
  beta03<-exp(summary(mod03)$coefficients[1,1])
  ic95a03<-exp(summary(mod03)$coefficients[1,1]-(z*summary(mod03)$coefficients[1,3]))
  ic95b03<-exp(summary(mod03)$coefficients[1,1]+(z*summary(mod03)$coefficients[1,3]))
  coef03<-ic_guapa2(guapa(beta03),guapa(ic95a03),guapa(ic95b03))
  pval03<-pval_guapa(summary(mod03)$coefficients[1,5])
  
  mod03_nl<-clogit(enfneuro~bs(aaa,degree=3)+strata(match)
                   +as.factor(grup_int)+as.factor(escolar)
                   +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                   +as.factor(apoe4),
                   data=datx)
  p_nonlin<-pval_guapa(lrtest(mod03,mod03_nl)[2,5])
  
  aaa<-datx[,vars02[i]]
  mod02<-clogit(enfneuro~datx[,vars01[i]]+strata(match)
                +as.factor(grup_int)+as.factor(escolar)
                +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                +as.factor(apoe4),
                data=datx)
  mod02_nl<-clogit(enfneuro~bs(aaa,degree=3)+strata(match)
                   +as.factor(grup_int)+as.factor(escolar)
                   +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)+imc
                   +as.factor(apoe4),
                   data=datx)
  
  ptemp<-termplot(mod02_nl,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,median(aaa,na.rm=TRUE))[1]
  center<-with(temp, y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./splines/",vars02[i],".jpg",sep="")
  labely<-c("Neurodegenerative disease risk (odds ratio)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  plot.data<-subset2(plot.data,"plot.data$lci>=0.1 & plot.data$uci<=10")
  
  p_lin2<-pval_guapa(summary(mod02)$coefficients[1,5])
  p_nonlin2<-pval_guapa(lrtest(mod02,mod02_nl)[2,5])
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,
             "\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(data=plot.data, aes_string(x=plot.data$x, y=plot.data$y)) + 
    geom_ribbon(aes_string(ymin=plot.data$lci, ymax=plot.data$uci), alpha=0.25, fill="Black") +
    geom_line(aes_string(x=plot.data$x, y=plot.data$y), color='Black') +
    geom_hline(yintercept=1, linetype=2) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    labs(x=vars03[i],y=labely) +
    annotate("text", x=max(plot.data$x,na.rm=TRUE)*0.98, y=max(plot.data$uci,na.rm=TRUE), label=leg, vjust=1, hjust=1, size=6) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  
  ggsave(filename=name, dpi=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03,p_nonlin))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02,
                         beta03,ic95a03,ic95b03))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (mod2)","pval (mod2)",
                 "OR (mod3)","pval (mod3)","pval-nonlin (adj)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./nested_case_control_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./nested_case_control_forestplots.csv",sep=";",col.names=NA)


### GENERATION OF FOREST PLOT ###

dat<-read.csv2("./nested_case_control_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)]),
                         as.matrix(dat[c(1,8,9,10)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,1,
             2,2,2,2,2,
             3,3,3,3,3)
dat$model<-with(dat,ifelse(model==1,"03. Not adjusted",
                           ifelse(model==2,"02. Adjusted",
                                  ifelse(model==3,"01. +imc",NA))))
dat$bmk<-c("01.\nLeptin","02.\nGhrelin","03.\nResistin","04.\nAdipsin",
           "05.\nFibroblast growth\nfactor-21",
           "01.\nLeptin","02.\nGhrelin","03.\nResistin","04.\nAdipsin",
           "05.\nFibroblast growth\nfactor-21",
           "01.\nLeptin","02.\nGhrelin","03.\nResistin","04.\nAdipsin",
           "05.\nFibroblast growth\nfactor-21")
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./nested_case_control_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(2)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./nested_case_control_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5.5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + 
    ylab("Risk of dementia\nper 1 SD higher levels of biomarkers at baseline\n(odds ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.25,0.5,1,2,4)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=13, scales="free_y") +
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
          axis.title.y=element_text(size=20),
          strip.text=element_text(size=8.5))
  
  ggsave("./nested_case_control_forestplot.tiff", units="px", width=9000, height=7500, dpi=1200, bg="transparent")
  figure
}


### SURVIVAL ANALYSES ###
#########################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_total.RData")

# IMPUTATION OF escolar AND apoe4 NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-with(dat_mice,ifelse(apoe4==9,NA,apoe4))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
round(sum(is.na(dat_mice$escolar))/dim(dat_mice)[1]*100,1) #1.8% NAs in escolar
round(sum(is.na(dat_mice$apoe4))/dim(dat_mice)[1]*100,1) #3% NAs in escolar

dat_mice[] <- lapply(dat_mice, function(x) if (is.factor(x) & !identical(x, dat_mice$escolar) & !identical(x, dat_mice$apoe4)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
dat_mice<-dat_mice[,c("id","escolar","apoe4")]
dat$escolar<-NULL
dat$apoe4<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)


dat$adjx<-1

z<-qnorm(1-0.05/2)

vars01<-c("z_imc","z_cintura")
vars02<-c("adjx","adjx")
vars03<-c("imc","cintura")
vars04<-c("Body mass index, kg/m2","Waist circumference, cm")

tab<-NULL
tab2<-NULL


for(i in 1:length(vars01))
  
{ 
  participants<-length(which(!is.na(dat[,vars01[i]])))
  dat2<-dat[!is.na(dat[,vars01[i]]),]
  aaa<-dat2[,vars01[i]]
  bbb<-dat2[,vars02[i]]
  
  mod01<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2),
               na.action=na.exclude, method="breslow", data=dat2)
  beta01<-intervals(mod01)[1,1]
  ic95a01<-intervals(mod01)[1,2]
  ic95b01<-intervals(mod01)[1,3]
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(intervals(mod01)[1,4])
  
  mod02<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2)
               +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
               +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)
               +bbb+as.factor(apoe4),
               na.action=na.exclude, method="breslow", data=dat2)
  beta02<-intervals(mod02)[1,1]
  ic95a02<-intervals(mod02)[1,2]
  ic95b02<-intervals(mod02)[1,3]
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(intervals(mod02)[1,4])
  
  aaa<-dat2[,vars03[i]]
  mod02<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2)
               +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
               +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)
               +as.factor(apoe4),
               na.action=na.exclude, method="breslow", data=dat2)
  mod02_nl<-coxph(Surv(toenfneuro, enfneuro)~pspline(aaa, df=4)+cluster(idcluster2)
                  +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
                  +as.factor(diabetes)+as.factor(hipercol)+as.factor(hta)+as.factor(tabaco)
                  +as.factor(apoe4),
                  na.action=na.exclude, method="breslow",data=dat2)
  p_nonlin<-pval_guapa(lrtest(mod02,mod02_nl)[2,5])
  
  ptemp<-termplot(mod02_nl,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,median(aaa,na.rm=TRUE))[1]
  center<-with(temp, y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./splines/",vars03[i],".jpg",sep="")
  labely<-c("Neurodegenerative disease risk (hazard ratio)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  plot.data<-subset2(plot.data,"plot.data$lci>=0.1 & plot.data$uci<=10")
  
  p_lin2<-pval_guapa(intervals(mod02)[1,4])
  p_nonlin2<-pval_guapa(lrtest(mod02,mod02_nl)[2,5])
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,
             "\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(data=plot.data, aes_string(x=plot.data$x, y=plot.data$y)) + 
    geom_ribbon(aes_string(ymin=plot.data$lci, ymax=plot.data$uci), alpha=0.25, fill="Black") +
    geom_line(aes_string(x=plot.data$x, y=plot.data$y), color='Black') +
    geom_hline(yintercept=1, linetype=2) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    labs(x=vars04[i],y=labely) +
    annotate("text", x=max(plot.data$x,na.rm=TRUE)*0.98, y=max(plot.data$uci,na.rm=TRUE), label=leg, vjust=1, hjust=1, size=6) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  
  ggsave(filename=name, dpi=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,p_nonlin))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (adj)","pval (adj)","non-linearity (pval)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./adiposity_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./adiposity_forestplots.csv",sep=";",col.names=NA)


### GENERATION OF FOREST PLOT (ONLY ADIPOSITY MEASURES) ###

dat<-read.csv2("./adiposity_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,
             2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-c("01.\nBody mass\nindex (+1 SD)","02.\nWaist circ.\n(+1 SD)",
           "01.\nBody mass\nindex (+1 SD)","02.\nWaist circ.\n(+1 SD)")
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./adiposity_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(1.5)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./adiposity_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5.5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + 
    ylab("Risk of dementia\n(hazard ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.25,0.5,1,2)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=13, scales="free_y") +
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
          axis.title.y=element_text(size=20),
          strip.text=element_text(size=6.5))
  
  ggsave("./adiposity_forestplot.tiff", units="px", width=9000, height=3000, dpi=1200, bg="transparent")
  figure
}



###############################
### ANALYSES: LIPID PROFILE ###
###############################

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Outputs/Lipids")
dir.create("./splines")


##############################################################################
### BIOMARKERS: LINYING MODIFICATION OF SURVIVAL ANALYSES FOR CASE-COHORTS ###
##############################################################################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis01_cc.RData")
dat$nodo<-with(dat,ifelse(nodo==12,11,nodo))

# IMPUTATION OF LIPID PROFILE NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hipertg","hta","apoe4",
                 "tabaco","imc","hdlc","ldlc","tg")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
round(sum(is.na(dat_mice$escolar))/dim(dat_mice)[1]*100,1) #1.0%
round(sum(is.na(dat_mice$ldlc))/dim(dat_mice)[1]*100,1) #4.2%
round(sum(is.na(dat_mice$hdlc))/dim(dat_mice)[1]*100,1) #2.2%
round(sum(is.na(dat_mice$tg))/dim(dat_mice)[1]*100,1) #1.6%


dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice<-dat_mice[,c("id","escolar","hdlc","ldlc","tg")]
dat$escolar<-NULL
dat$z_hdlc<-as.numeric(with(dat,scale(dat$hdlc)))
dat$z_ldlc<-as.numeric(with(dat,scale(dat$ldlc)))
dat$z_logtg<-as.numeric(with(dat,scale(log(dat$tg))))
dat$hdlc<-NULL
dat$ldlc<-NULL
dat$tg<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


# ANALYSES #

vars01<-c("z_eflnor01","z_vcam01","z_h2401","z_apoe01")

ntotal<-7447
z<-qnorm(1-0.05/2)
vars<-c("id","toenfneuro","enfneuro","subcohort","edad","sexo","grup_int","escolar",
        "diabetes","hipercol","hipertg","hta","tabaco","imc","apoe4","hdlc","ldlc","tg")

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
             +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)
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
  
  mod03<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]]
             +edad+as.factor(sexo)+as.factor(grup_int)+as.factor(escolar)
             +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)
             +as.factor(hta)+as.factor(tabaco)+imc+as.factor(apoe4)+hdlc+ldlc+tg,
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta03<-exp(summary(mod03)$coefficients[1,1])
  ic95a03<-exp(summary(mod03)$coefficients[1,1]-(z*summary(mod03)$coefficients[1,2]))
  ic95b03<-exp(summary(mod03)$coefficients[1,1]+(z*summary(mod03)$coefficients[1,2]))
  coef03<-ic_guapa2(guapa(beta03),guapa(ic95a03),guapa(ic95b03))
  pval03<-pval_guapa(summary(mod03)$coefficients[1,4])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02,beta03,ic95a03,ic95b03))
}

colnames(tab)<-c("HR (raw)","pval (raw)","HR (mod2)","pval (mod2)",
                 "HR (mod3)","pval (mod3)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./case_cohort_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./case_cohort_forestplots.csv",sep=";",col.names=NA)


dat<-read.csv2("./case_cohort_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)]),
                         as.matrix(dat[c(1,8,9,10)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,
             2,2,2,2,
             3,3,3,3)
dat$model<-with(dat,ifelse(model==1,"03. Not adjusted",
                           ifelse(model==2,"02. Adjusted",
                                  ifelse(model==3,"01. +lipids",NA))))
dat$bmk<-with(dat,ifelse(X=="z_eflnor01","01.\nEfflux\n(SH-SY5Y cells)",
                         ifelse(X=="z_vcam01","021.\nVCAM1 release\n(HBEC-5i cells)",
                                ifelse(X=="z_h2401","03.\n24S-hydroxy-\ncholesterol",
                                       ifelse(X=="z_apoe01","04.\nApolipoprotein\nE",NA)))))
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
          strip.text=element_text(size=6.5))
  
  namefile<-paste("./case_cohort_",vars01[i],".tiff",sep="")
  ggsave(filename=namefile, units="px", width=9000, height=6000, dpi=1200, bg="transparent")
  figure
}


### CAMBIOS A 1 ANY ###

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis01_cc.RData")
dat$nodo<-with(dat,ifelse(nodo==12,11,nodo))

# IMPUTATION OF LIPID PROFILE NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hipertg","hta","apoe4",
                 "tabaco","imc","hdlc","ldlc","tg")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))

dat_mice[] <- lapply(dat_mice, function(x) if(is.factor(x) & !identical(x, dat_mice$escolar)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice<-dat_mice[,c("id","escolar","hdlc","ldlc","tg")]
dat$escolar<-NULL
dat$hdlc<-NULL
dat$ldlc<-NULL
dat$tg<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)
dat$escolar<-with(dat,ifelse(escolar==3,2,escolar))


# ANALYSES #

vars01<-c("z_eflnor01","z_vcam01","z_h2401","z_apoe01")
vars02<-c("z_eflnor03","z_vcam03","z_h2403","z_apoe03")

ntotal<-7447
z<-qnorm(1-0.05/2)
vars<-c("id","toenfneuro","enfneuro","subcohort","edad","sexo","grup_int","escolar",
        "diabetes","hipercol","hipertg","hta","tabaco","imc","apoe4","hdlc","ldlc","tg")

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
             +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)
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
  
  mod03<-cch(Surv(toenfneuro,enfneuro) ~datx[,vars01[i]]+datx[,vars02[i]]
             +edad+as.factor(sexo)+as.factor(grup_int)+as.factor(escolar)
             +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)
             +as.factor(hta)+as.factor(tabaco)+imc+as.factor(apoe4)+hdlc+ldlc+tg,
             subcoh= ~subcohort,
             data=datx,
             id= ~id,
             cohort.size=ntotal,
             method="LinYing")
  beta03<-exp(summary(mod03)$coefficients[1,1])
  ic95a03<-exp(summary(mod03)$coefficients[1,1]-(z*summary(mod03)$coefficients[1,2]))
  ic95b03<-exp(summary(mod03)$coefficients[1,1]+(z*summary(mod03)$coefficients[1,2]))
  coef03<-ic_guapa2(guapa(beta03),guapa(ic95a03),guapa(ic95b03))
  pval03<-pval_guapa(summary(mod03)$coefficients[1,4])
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02,beta03,ic95a03,ic95b03))
}

colnames(tab)<-c("HR (raw)","pval (raw)","HR (mod2)","pval (mod2)",
                 "HR (mod3)","pval (mod3)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./case_cohort_1year_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./case_cohort_1year_forestplots.csv",sep=";",col.names=NA)


dat<-read.csv2("./case_cohort_1year_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)]),
                         as.matrix(dat[c(1,8,9,10)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,
             2,2,2,2,
             3,3,3,3)
dat$model<-with(dat,ifelse(model==1,"03. Not adjusted",
                           ifelse(model==2,"02. Adjusted",
                                  ifelse(model==3,"01. +lipids",NA))))
dat$bmk<-with(dat,ifelse(X=="z_eflnor01","01.\nEfflux\n(SH-SY5Y cells)",
                         ifelse(X=="z_vcam01","021.\nVCAM1 release\n(HBEC-5i cells)",
                                ifelse(X=="z_h2401","03.\n24S-hydroxy-\ncholesterol",
                                       ifelse(X=="z_apoe01","04.\nApolipoprotein\nE",NA)))))
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
          strip.text=element_text(size=6.5))
  
  namefile<-paste("./case_cohort_1year_",vars01[i],".tiff",sep="")
  ggsave(filename=namefile, units="px", width=9000, height=6000, dpi=1200, bg="transparent")
  figure
}


### SURVIVAL ANALYSES ###
#########################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_total.RData")

# IMPUTATION OF escolar AND apoe4 NAs #

dat_mice<-dat[,c("id","enfneuro","sexo","edad","escolar","grup_int","diabetes","hipercol","hta","apoe4",
                 "tabaco","imc")]
dat_mice$escolar<-with(dat_mice,ifelse(escolar==9,NA,escolar))
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-with(dat_mice,ifelse(apoe4==9,NA,apoe4))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
round(sum(is.na(dat_mice$escolar))/dim(dat_mice)[1]*100,1) #1.8% NAs in escolar
round(sum(is.na(dat_mice$apoe4))/dim(dat_mice)[1]*100,1) #3% NAs in escolar

dat_mice[] <- lapply(dat_mice, function(x) if (is.factor(x) & !identical(x, dat_mice$escolar) & !identical(x, dat_mice$apoe4)) as.character(x) else x)
imputed_data<-missForest(dat_mice)
dat_mice<-imputed_data$ximp
dat_mice$escolar<-factor(dat_mice$escolar,levels=c(1,2,3))
dat_mice$apoe4<-factor(dat_mice$apoe4,levels=c(0,1))
dat_mice<-dat_mice[,c("id","escolar","apoe4")]
dat$escolar<-NULL
dat$apoe4<-NULL
dat<-merge2(dat,dat_mice,by.id=c("id"),all.x=TRUE,sort=FALSE)


dat<-dat[!is.na(dat$hdlc) | !is.na(dat$ldlc) | !is.na(dat$tg),]

z<-qnorm(1-0.05/2)

vars01<-c("z_tc","z_ldlc","z_hdlc","z_logtg")
vars02<-c("tc","ldlc","hdlc","tg")
vars04<-c("Total cholesterol, mg/dL","LDL cholesterol, mg/dL","HDL cholesterol, mg/dL","Triglycerides, mg/dL")

tab<-NULL
tab2<-NULL


for(i in 1:length(vars01))
  
{ 
  participants<-length(which(!is.na(dat[,vars01[i]])))
  dat2<-dat[!is.na(dat[,vars01[i]]),]
  aaa<-dat2[,vars01[i]]

  mod01<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2),
               na.action=na.exclude, method="breslow", data=dat2)
  beta01<-intervals(mod01)[1,1]
  ic95a01<-intervals(mod01)[1,2]
  ic95b01<-intervals(mod01)[1,3]
  coef01<-ic_guapa2(guapa(beta01),guapa(ic95a01),guapa(ic95b01))
  pval01<-pval_guapa(intervals(mod01)[1,4])
  
  mod02<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2)
               +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
               +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)+as.factor(hta)+as.factor(tabaco)
               +as.factor(apoe4)+imc,
               na.action=na.exclude, method="breslow", data=dat2)
  beta02<-intervals(mod02)[1,1]
  ic95a02<-intervals(mod02)[1,2]
  ic95b02<-intervals(mod02)[1,3]
  coef02<-ic_guapa2(guapa(beta02),guapa(ic95a02),guapa(ic95b02))
  pval02<-pval_guapa(intervals(mod02)[1,4])
  
  aaa<-dat2[,vars02[i]]
  mod02<-coxph(Surv(toenfneuro, enfneuro)~aaa+cluster(idcluster2)
               +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
               +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)+as.factor(hta)+as.factor(tabaco)
               +as.factor(apoe4)+imc,
               na.action=na.exclude, method="breslow", data=dat2)
  mod02_nl<-coxph(Surv(toenfneuro, enfneuro)~pspline(aaa, df=4)+cluster(idcluster2)
                  +as.factor(grup_int)+edad+as.factor(sexo)+as.factor(escolar)+as.factor(nodo)
                  +as.factor(diabetes)+as.factor(hipercol)+as.factor(hipertg)+as.factor(hta)+as.factor(tabaco)
                  +as.factor(apoe4)+imc,
                  na.action=na.exclude, method="breslow",data=dat2)
  p_nonlin<-pval_guapa(lrtest(mod02,mod02_nl)[2,5])
  
  ptemp<-termplot(mod02_nl,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,median(aaa,na.rm=TRUE))[1]
  center<-with(temp, y[x==value])
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./splines/",vars02[i],".jpg",sep="")
  labely<-c("Neurodegenerative disease risk (hazard ratio)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  plot.data<-subset2(plot.data,"plot.data$lci>=0.1 & plot.data$uci<=10")
  
  p_lin2<-pval_guapa(intervals(mod02)[1,4])
  p_nonlin2<-pval_guapa(lrtest(mod02,mod02_nl)[2,5])
  p_lin2<-ifelse(p_lin2=="<0.001"," < 0.001",
                 ifelse(p_lin2=="<0.00001"," < 0.00001",paste(" = ",p_lin2,sep="")))
  p_nonlin2<-ifelse(p_nonlin2=="<0.001"," < 0.001",
                    ifelse(p_nonlin2=="<0.00001"," < 0.00001",paste(" = ",p_nonlin2,sep="")))
  leg<-paste("p-value for linearity",p_lin2,
             "\np-value for non-linearity",p_nonlin2,sep="")
  
  figure<-ggplot(data=plot.data, aes_string(x=plot.data$x, y=plot.data$y)) + 
    geom_ribbon(aes_string(ymin=plot.data$lci, ymax=plot.data$uci), alpha=0.25, fill="Black") +
    geom_line(aes_string(x=plot.data$x, y=plot.data$y), color='Black') +
    geom_hline(yintercept=1, linetype=2) +
    theme_bw() +
    scale_x_continuous(expand=c(0,0)) +
    labs(x=vars04[i],y=labely) +
    annotate("text", x=max(plot.data$x,na.rm=TRUE)*0.98, y=max(plot.data$uci,na.rm=TRUE), label=leg, vjust=1, hjust=1, size=6) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank())  
  
  ggsave(filename=name, dpi=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,p_nonlin))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("OR (raw)","pval (raw)","OR (adj)","pval (adj)","non-linearity (pval)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./lipids_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./lipids_forestplots.csv",sep=";",col.names=NA)


### GENERATION OF FOREST PLOT ###

dat<-read.csv2("./lipids_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,
             2,2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-c("01.\nTotal cholesterol\n(+1 SD)","02.\nLDL cholesterol\n(+1 SD)",
           "03.\nHDL cholesterol\n(+1 SD)","04.\nTriglycerides\n(log) (+1 SD)",
           "01.\nTotal cholesterol\n(+1 SD)","02.\nLDL cholesterol\n(+1 SD)",
           "03.\nHDL cholesterol\n(+1 SD)","04.\nTriglycerides\n(log) (+1 SD)")
dat$coef<-ic_guapa2(guapa(as.numeric(dat$est)),guapa(as.numeric(dat$lo)),guapa(as.numeric(dat$hi)))
write.table(dat,file="./lipids_forestplot.csv",sep=";",col.names=TRUE,row.names=FALSE)

vars01<-c("forestplot")
vars02<-c(1.5)

for(i in 1:length(vars01))
  
{
  dat<-read.csv2("./lipids_forestplot.csv",header=TRUE,sep=";",dec=".")
  
  figure<-ggplot(data=dat,
                 aes(x=model, y=est, ymin=lo, ymax=hi)) +
    geom_hline(aes(fill=model), yintercept=1, linetype=2) +
    geom_pointrange(aes(col=model), size=0.7, shape=15) +
    geom_text(data=dat, size=5.5, aes(y=max(hi)*vars02[i], x=model, label=coef, hjust='inward')) +
    xlab(" ") + 
    ylab("Risk of dementia\n(hazard ratio, 95% confidence interval)") +
    scale_y_continuous(trans = log2_trans(), breaks=c(0.85,1,1.25)) +
    geom_errorbar(aes(ymin=lo, ymax=hi, col=model), width=0.6, cex=1) + 
    facet_wrap(~bmk, strip.position="left", nrow=13, scales="free_y") +
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
          axis.title.y=element_text(size=20),
          strip.text=element_text(size=6.5))
  
  ggsave("./lipids_forestplot.tiff", units="px", width=9000, height=5400, dpi=1200, bg="transparent")
  figure
}



##############################
### ANALYSES: INFLAMMATION ###
##############################

setwd("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Outputs/Inflammation")
dir.create("./splines")


##############################################################################
### BIOMARKERS: LINYING MODIFICATION OF SURVIVAL ANALYSES FOR CASE-COHORTS ###
##############################################################################

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis01_cc.RData")
dat$nodo<-with(dat,ifelse(nodo==12,11,nodo))

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

vars01<-c("z_log_ldh01","z_c301","z_antitrip01","z_a2m01","z_ppy01")
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
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02))
  tab2<-rbind(tab2,cbind(beta01,ic95a01,ic95b01,beta02,ic95a02,ic95b02))
}

colnames(tab)<-c("HR (raw)","pval (raw)","HR (mod2)","pval (mod2)")
rownames(tab)<-vars01
rownames(tab2)<-vars01
write.table(tab,file="./case_cohort_results.csv",sep=";",col.names=NA)
write.table(tab2,file="./case_cohort_forestplots.csv",sep=";",col.names=NA)


dat<-read.csv2("./case_cohort_forestplots.csv",header=TRUE,sep=";",dec=".")
dat<-as.data.frame(rbind(as.matrix(dat[1:4]),as.matrix(dat[c(1,5,6,7)])))
names(dat)<-c("X","est","lo","hi")
dat$model<-c(1,1,1,1,1,
             2,2,2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-with(dat,ifelse(X=="z_log_ldh01","01.\nLactate\ndehydrogenase",
                         ifelse(X=="z_c301","02.\nComplement\ncomponent 3",
                                ifelse(X=="z_antitrip01","03.\nAlpha-1\nantitripsin",
                                       ifelse(X=="z_a2m01","04.\nAlpha-2\nmacroglobulin",
                                              ifelse(X=="z_ppy01","05.\nPancreatic\nprohormone",NA))))))
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
  ggsave(filename=namefile, units="px", width=9000, height=5700, dpi=1200, bg="transparent")
  figure
}


### CAMBIOS A 1 ANY ###

load("C:/Users/ALHE/OneDrive - Folkehelseinstituttet/Articles/PREDIMED_dementia/Data/predimed_dem_fis01_cc.RData")
dat$nodo<-with(dat,ifelse(nodo==12,11,nodo))

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

vars01<-c("z_log_ldh01","z_c301","z_antitrip01","z_a2m01","z_ppy01")
vars02<-c("z_log_ldh03","z_c303","z_antitrip03","z_a2m03","z_ppy03")

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
dat$model<-c(1,1,1,1,1,
             2,2,2,2,2)
dat$model<-with(dat,ifelse(model==1,"02. Not adjusted",
                           ifelse(model==2,"01. Adjusted",NA)))
dat$bmk<-with(dat,ifelse(X=="z_log_ldh01","01.\nLactate\ndehydrogenase",
                         ifelse(X=="z_c301","02.\nComplement\ncomponent 3",
                                ifelse(X=="z_antitrip01","03.\nAlpha-1\nantitripsin",
                                       ifelse(X=="z_a2m01","04.\nAlpha-2\nmacroglobulin",
                                              ifelse(X=="z_ppy01","05.\nPancreatic\nprohormone",NA))))))
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
  ggsave(filename=namefile, units="px", width=9000, height=5700, dpi=1200, bg="transparent")
  figure
}

