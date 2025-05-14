library(survival)
library(survidm)
library(Rcpp)
library(RcppRoll)
library(ggplot2)
library(ggpubr)

ivus=function(times,VUS_result,cox12_marginal,cox13_marginal,cox23_marginal) {
  modvecs=get_modvecs(dat,cox12_marginal,cox13_marginal,cox23_marginal)
  result=P(times,modvecs$baselambda12,modvecs$explp12,
           modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)[[1]]
  indices=which(times %in% VUS_result[[4]])
  P1=result[1,indices]
  P2=result[2,indices]
  P3=result[3,indices]
  P3_temp=result[3,indices]-c(0,result[3,indices][-length(indices)])
  
  return(sum(VUS_result[[1]]*P1*P2*P3_temp,na.rm=T)/sum(P1*P2*P3_temp,na.rm=T))
}

explp_calc=function(df,times) {
  return(rep(df$explp[1],length(times)))
}

get_modvecs=function(dat,cox12,cox13,cox23) {
  times=sort(unique(dat$time[dat$event==1]))
  BaseLambda12=basehaz(cox12,centered=F)
  BaseLambda12_fun=stepfun(BaseLambda12$time,c(0,BaseLambda12$hazard))
  BaseLambda13=basehaz(cox13,centered=F)
  BaseLambda13_fun=stepfun(BaseLambda13$time,c(0,BaseLambda13$hazard))
  BaseLambda23=basehaz(cox23,centered=F)
  BaseLambda23_fun=stepfun(BaseLambda23$time,c(0,BaseLambda23$hazard))
  baselambda12=BaseLambda12_fun(times)
  baselambda13=BaseLambda13_fun(times)
  baselambda23=BaseLambda23_fun(times)
  explp12_df=data.frame(id=dat$id,
                        explp=predict(cox12,newdata=dat,type="risk"))
  explp13_df=data.frame(id=dat$id,
                        explp=predict(cox13,newdata=dat,type="risk"))
  explp23_df=data.frame(id=dat$id,
                        explp=predict(cox23,newdata=dat,type="risk"))
  
  explp12=do.call(rbind,lapply(split(explp12_df,explp12_df$id),explp_calc,times=times))
  explp13=do.call(rbind,lapply(split(explp13_df,explp13_df$id),explp_calc,times=times))
  explp23=do.call(rbind,lapply(split(explp23_df,explp23_df$id),explp_calc,times=times))
  
  return(list(baselambda12=baselambda12,baselambda13=baselambda13,baselambda23=baselambda23,
              explp12=explp12,explp13=explp13,explp23=explp23))
}

smooth=function(HUM_ID,window=50,method="ID") {
  hums=HUM_ID[[1]]
  conc=HUM_ID[[2]]
  comp=HUM_ID[[3]]
  time=HUM_ID[[4]]
  hums=hums[order(time)]
  conc=conc[order(time)]
  comp=comp[order(time)]
  time=time[order(time)]
  vals=roll_sum(conc,window)/roll_sum(comp,window)
  vals_time=time[(window/2):(length(time)-(window/2))]
  vals=roll_sum(conc,window)/roll_sum(comp,window)
  df_ID=data.frame(vals_time,vals,method=rep(method,length(vals)))
  return(df_ID)
}

colonIDM=colonIDM[colonIDM$rx=="Obs",]
colonIDM$time=colonIDM$Stime
colonIDM_trans23=colonIDM[colonIDM$event1==1,]
dat=colonIDM
dat$id=1:dim(dat)[1]
dat_trans23=colonIDM_trans23
ill_times=dat$time1
ill_times[dat$time1==dat$time & dat$event==1]=Inf
death_times=dat$time
death_times[dat$event==0]=Inf
censoring_times=dat$time
censoring_times[dat$event==1]=Inf
G=rep(0,length(death_times))
for(i in 1:length(death_times)) {
  if(death_times[i]!=Inf) {
    G[i]=summary(survfit(Surv(dat$time,1-dat$event)~1),times=death_times[i])[[6]]
  }
}

#Marginal Models
cox12_marginal=coxph(Surv(time1,event1)~1,data=dat)
cox13_marginal=coxph(Surv(time1,event)~1,data=dat)
cox23_marginal=coxph(Surv(time1,time,event)~1,data=dat_trans23)

#Model A
cox12=coxph(Surv(time1,event1)~age+sex,data=dat)
cox13=coxph(Surv(time1,event)~age+sex,data=dat)
cox23=coxph(Surv(time1,time,event)~age+sex,data=dat_trans23)
modvecs=get_modvecs(dat,cox12,cox13,cox23)
HUM_ID_A=VUS_ID(sort(unique(dat$time[dat$event==1])),ill_times,death_times,censoring_times,modvecs$baselambda12,modvecs$explp12,
                modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
df_ID_A=smooth(HUM_ID_A)
ivusA=ivus(sort(unique(dat$time[dat$event==1])),HUM_ID_A,cox12_marginal,cox13_marginal,cox23_marginal)
cstatsA=c(concordance(cox12)$concordance,concordance(cox13)$concordance,concordance(cox23)$concordance)
HUM_CD_A=VUS_CD(sort(unique(dat$time[dat$event==1])),ill_times,death_times,censoring_times,G,modvecs$baselambda12,modvecs$explp12,
                modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
df_CD_A=smooth(HUM_CD_A,method="CD")

#Model B
cox12=coxph(Surv(time1,event1)~age+sex+obstruct+adhere+nodes+differ+extent+surg,data=dat)
cox13=coxph(Surv(time1,event)~age+sex+obstruct+adhere+nodes+differ+extent+surg,data=dat)
cox23=coxph(Surv(time1,time,event)~age+sex+obstruct+adhere+nodes+differ+extent+surg,data=dat_trans23)
modvecs=get_modvecs(dat,cox12,cox13,cox23)
HUM_ID_B=VUS_ID(sort(unique(dat$time[dat$event==1])),ill_times,death_times,censoring_times,modvecs$baselambda12,modvecs$explp12,
                modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
df_ID_B=smooth(HUM_ID_B)
ivusB=ivus(sort(unique(dat$time[dat$event==1])),HUM_ID_B,cox12_marginal,cox13_marginal,cox23_marginal)
cstatsB=c(concordance(cox12)$concordance,concordance(cox13)$concordance,concordance(cox23)$concordance)
HUM_CD_B=VUS_CD(sort(unique(dat$time[dat$event==1])),ill_times,death_times,censoring_times,G,modvecs$baselambda12,modvecs$explp12,
                modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
df_CD_B=smooth(HUM_CD_B,method="CD")

#Results
cstats=data.frame(cstatsA,cstatsB)
ivus_vals=data.frame(ivusA,ivusB)
cstats
ivus_vals

write.csv(cstats,"//Users//nicholashartman//Documents//VUS//cstats.csv")
write.csv(ivus_vals,"//Users//nicholashartman//Documents//VUS//ivus_vals.csv")
save(HUM_ID_A,file="//Users//nicholashartman//Documents//VUS//HUM_ID_A.Rdata")
save(HUM_CD_A,file="//Users//nicholashartman//Documents//VUS//HUM_CD_A.Rdata")
save(HUM_ID_B,file="//Users//nicholashartman//Documents//VUS//HUM_ID_B.Rdata")
save(HUM_CD_B,file="//Users//nicholashartman//Documents//VUS//HUM_CD_B.Rdata")
write.csv(rbind(df_ID_A,df_CD_A),"//Users//nicholashartman//Documents//VUS//df_A.csv",row.names=F)
write.csv(rbind(df_ID_B,df_CD_B),"//Users//nicholashartman//Documents//VUS//df_B.csv",row.names=F)
df_ID_A$model="Model A"
df_ID_B$model="Model B"
df_CD_A$model="Model A"
df_CD_B$model="Model B"
write.csv(rbind(df_ID_A,df_ID_B),"//Users//nicholashartman//Documents//VUS//df_ID.csv",row.names=F)
write.csv(rbind(df_CD_A,df_CD_B),"//Users//nicholashartman//Documents//VUS//df_CD.csv",row.names=F)

df=read.csv("//Users//nicholashartman//Documents//VUS//df_ID.csv")
g=ggplot(df)+geom_line(aes(x=vals_time,y=vals,colour=model,linetype=model),lwd=1.2)
g=g+scale_colour_manual(values=c("black","dodgerblue3"))+scale_linetype_manual(values=c("dotted","dashed"))
g=g+geom_hline(aes(yintercept=1/6),colour="red",lwd=1.2)
g=g+theme_classic()+ylim(0,0.4)+xlab("Time (days)")+ylab("")
g=g+theme(text=element_text(size=24),legend.position="top",legend.key.width=unit(2.3,"cm"),legend.title=element_blank())
g=g+ggtitle("ID")

df=read.csv("//Users//nicholashartman//Documents//VUS//df_CD.csv")
g2=ggplot(df)+geom_line(aes(x=vals_time,y=vals,colour=model,linetype=model),lwd=1.2)
g2=g2+scale_colour_manual(values=c("black","dodgerblue3"))+scale_linetype_manual(values=c("dotted","dashed"))
g2=g2+geom_hline(aes(yintercept=1/6),colour="red",lwd=1.2)
g2=g2+theme_classic()+ylim(0,0.4)+xlab("Time (days)")+ylab("HUM")
g2=g2+theme(text=element_text(size=24),legend.position="top",legend.key.width=unit(2.3,"cm"),legend.title=element_blank())
g2=g2+ggtitle("CD")

ggarrange(g2,g,nrow=1,common.legend=T)