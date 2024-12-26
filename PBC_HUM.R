library(survival)
library(Rcpp)
library(RcppRoll)
library(ggplot2)
library(ggpubr)
sourceCpp("PBC_Functions.cpp")

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

################
#Data Processing
################
dat=pbc[pbc$edema==0,]

pbc2=tmerge(pbc, pbc, id=id, death=event(time, status)) #set range
pbc2=tmerge(pbc2, pbcseq, id=id, ascites=tdc(day, ascites),
            bili=tdc(day, bili), albumin=tdc(day, albumin),
            protime=tdc(day, protime), alk.phos=tdc(day, alk.phos),
            edema=tdc(day,edema))
ids=pbc2$id[pbc2$edema==0 & pbc2$tstart==0]
pbc2=pbc2[pbc2$id %in% ids,]
dat=dat[dat$id %in% ids,]

get_edema_times=function(df) {
  out=ifelse(sum(df$edema) > 0,min(df$tstart[df$edema > 0]),max(df$tstop))
  return(out)
}
get_edema_stats=function(df) {
  out=ifelse(sum(df$edema) > 0,1,0)
  return(out)
}
dat$edema_time=sapply(split(pbc2,pbc2$id),get_edema_times)
dat$edema_stat=sapply(split(pbc2,pbc2$id),get_edema_stats)

##########################
#Prepare Modeling Datasets
##########################
dat$death=as.numeric(dat$status > 0)
dat$death1=as.numeric(dat$edema_stat==0 & dat$death==1)
dat_trans23=dat[dat$edema_stat==1,]

#################################
#Functions to Prepare Risk Scores
#################################
explp_calc=function(df,times) {
  return(rep(df$explp[1],length(times)))
}

get_modvecs=function(dat,cox12,cox13,cox23) {
  times=sort(unique(dat$time[dat$death==1]))
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

#############
#Estimate HUM
#############

#Marginal Models
cox12_marginal=coxph(Surv(edema_time,edema_stat)~1,data=dat)
cox13_marginal=coxph(Surv(edema_time,death1)~1,data=dat)
cox23_marginal=coxph(Surv(time,death)~1,data=dat_trans23)

#Model A
cox12=coxph(Surv(edema_time,edema_stat)~age+factor(sex)+stage,data=dat)
cox13=coxph(Surv(edema_time,death1)~age+factor(sex)+stage,data=dat)
cox23=coxph(Surv(time,death)~age+factor(sex)+stage,data=dat_trans23)

modvecs=get_modvecs(dat,cox12,cox13,cox23)
ill_times=dat$edema_time
ill_times[dat$edema_time==dat$time & dat$edema_stat==0]=Inf
death_times=dat$time
death_times[dat$death==0]=Inf
censoring_times=dat$time
censoring_times[dat$death==1]=Inf
G=rep(0,length(death_times))
for(i in 1:length(death_times)) {
  if(death_times[i]!=Inf) {
    G[i]=summary(survfit(Surv(dat$time,1-dat$death)~1),times=death_times[i])[[6]]
  }
}
HUM_ID=VUS_ID(sort(unique(dat$time[dat$death==1])),ill_times,death_times,censoring_times,modvecs$baselambda12,modvecs$explp12,
              modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
ivusA=ivus(sort(unique(dat$time[dat$death==1])),HUM_ID,cox12_marginal,cox13_marginal,cox23_marginal)
cstatsA=c(concordance(cox12)$concordance,concordance(cox13)$concordance,concordance(cox23)$concordance)
HUM_CD=VUS_CD(sort(unique(dat$time[dat$death==1])),ill_times,death_times,censoring_times,G,modvecs$baselambda12,modvecs$explp12,
              modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)

window=10
hums=HUM_ID[[1]]
conc=HUM_ID[[2]]
comp=HUM_ID[[3]]
time=HUM_ID[[4]]
hums=hums[order(time)]
conc=conc[order(time)]
comp=comp[order(time)]
time=time[order(time)]
vals_time=time[5:(length(time)-5)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
df_ID=data.frame(vals_time,vals,method=rep("ID",length(vals)))

hums=HUM_CD[[1]]
conc=HUM_CD[[2]]
comp=HUM_CD[[3]]
time=HUM_CD[[4]]
hums=hums[order(time)]
conc=conc[order(time)]
comp=comp[order(time)]
time=time[order(time)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
vals_time=time[5:(length(time)-5)]
df_CD=data.frame(vals_time,vals,method=rep("CD",length(vals)))

df=rbind(df_ID,df_CD)
write.csv(df,"pbc_hum_mod1.csv")

df=read.csv("pbc_hum_mod1.csv")
g=ggplot(df)+geom_line(aes(x=vals_time,y=vals,colour=method,linetype=method),lwd=1.2)
g=g+scale_colour_manual(values=c("black","dodgerblue3"))+scale_linetype_manual(values=c("dotted","dashed"))
g=g+geom_hline(aes(yintercept=1/6),colour="red",lwd=1.2)
g=g+theme_classic()+ylim(0,0.6)+xlab("Time (days)")+ylab("HUM")
g=g+theme(text=element_text(size=24),legend.position="top",legend.key.width=unit(2.3,"cm"),legend.title=element_blank())
g=g+ggtitle("Model A")

#Model B 
cox12=coxph(Surv(edema_time,edema_stat)~age+factor(sex)+stage+bili+albumin+protime,data=dat)
cox13=coxph(Surv(edema_time,death1)~age+factor(sex)+stage+bili+albumin+protime,data=dat)
cox23=coxph(Surv(time,death)~age+factor(sex)+stage+bili+albumin+protime,data=dat_trans23)

modvecs=get_modvecs(dat,cox12,cox13,cox23)
HUM_ID=VUS_ID(sort(unique(dat$time[dat$death==1])),ill_times,death_times,censoring_times,modvecs$baselambda12,modvecs$explp12,
              modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)
ivusB=ivus(sort(unique(dat$time[dat$death==1])),HUM_ID,cox12_marginal,cox13_marginal,cox23_marginal)
cstatsB=c(concordance(cox12)$concordance,concordance(cox13)$concordance,concordance(cox23)$concordance)
HUM_CD=VUS_CD(sort(unique(dat$time[dat$death==1])),ill_times,death_times,censoring_times,G,modvecs$baselambda12,modvecs$explp12,
              modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)

window=10
hums=HUM_ID[[1]]
conc=HUM_ID[[2]]
comp=HUM_ID[[3]]
time=HUM_ID[[4]]
hums=hums[order(time)]
conc=conc[order(time)]
comp=comp[order(time)]
time=time[order(time)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
vals_time=time[5:(length(time)-5)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
df_ID=data.frame(vals_time,vals,method=rep("ID",length(vals)))

hums=HUM_CD[[1]]
conc=HUM_CD[[2]]
comp=HUM_CD[[3]]
time=HUM_CD[[4]]
hums=hums[order(time)]
conc=conc[order(time)]
comp=comp[order(time)]
time=time[order(time)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
vals_time=time[5:(length(time)-5)]
vals=roll_sum(conc,window)/roll_sum(comp,window)
df_CD=data.frame(vals_time,vals,method=rep("CD",length(vals)))

df=rbind(df_ID,df_CD)
write.csv(df,"pbc_hum_mod2.csv")

df=read.csv("pbc_hum_mod2.csv")
g2=ggplot(df)+geom_line(aes(x=vals_time,y=vals,colour=method,linetype=method),lwd=1.2)
g2=g2+scale_colour_manual(values=c("black","dodgerblue3"))+scale_linetype_manual(values=c("dotted","dashed"))
g2=g2+geom_hline(aes(yintercept=1/6),colour="red",lwd=1.2)
g2=g2+theme_classic()+ylim(0,0.6)+xlab("Time (days)")+ylab("")
g2=g2+theme(text=element_text(size=24),legend.position="top",legend.key.width=unit(2.3,"cm"),legend.title=element_blank())
g2=g2+ggtitle("Model B")

ggarrange(g,g2,nrow=1,common.legend=T)
ggsave("PBC.pdf",width=15,height=7)

cstats=data.frame(cstatsA,cstatsB)
ivus_vals=data.frame(ivusA,ivusB)
write.csv(cstats,"cstats.csv",row.names=F)
write.csv(ivus_vals,"ivus_vals.csv",row.names=F)




