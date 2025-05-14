library(survival)
library(simsurv)
library(Rcpp)
sourceCpp("Functions.cpp")

get_data=function(n,tmax) {
  Z1=rnorm(n)
  Z2=rnorm(n)
  Z=data.frame(Z1,Z2)
  betas=c(0.5,0.5)
  names(betas)=c("Z1","Z2")
  tde=c(0,-0.5)
  names(tde)=c("Z1","Z2")
  ill=simsurv(dist="exponential",lambdas=0.025,x=Z,betas=betas,tde=tde,tdefunction=function(x){as.numeric(x > tmax/2)},maxt=tmax)
  names(ill)=c("id","ill_time","ill")
  betas=c(-0.5,-0.5)
  names(betas)=c("Z1","Z2")
  tde=c(0,0.5)
  names(tde)=c("Z1","Z2")
  death=simsurv(dist="exponential",lambdas=0.025,x=Z,betas=betas,tde=tde,tdefunction=function(x){as.numeric(x > tmax/2)},maxt=tmax)
  names(death)=c("id","death_time","death")
  betas=c(-0.5,-0.5)
  names(betas)=c("Z1","Z2")
  tde=c(0,0.5)
  names(tde)=c("Z1","Z2")
  post_death=simsurv(dist="exponential",lambdas=0.025,x=Z,betas=betas,tde=tde,tdefunction=function(x){as.numeric(x > tmax/2)},maxt=tmax)
  names(post_death)=c("id","post_death_time","post_death")
  result=merge(merge(ill,death,by="id"),post_death,by="id")
  result$ill=as.numeric(result$ill_time < result$death_time)
  result$ill_time=pmin(result$ill_time,result$death_time)
  result$death_time[result$ill==1]=result$ill_time[result$ill==1]+result$post_death_time[result$ill==1]
  result$death1=as.numeric(result$ill==0 & result$death==1)
  result=cbind(result,Z)
  return(result)
}

explp_calc=function(df,cox,times,ttfun,type) {
  if(type=="A") {
    explp=rep(predict(cox,newdata=df,type="risk")[1],length(times))
  } else {
    coef_cox=coef(cox)
    explp=exp(coef_cox[1]*df$Z1[1]+coef_cox[2]*df$Z2[1]+coef_cox[3]*df$Z2[1]*ttfun(times))
  }
  return(explp)
}

get_modvecs=function(dat,cox12,cox13,cox23,ttfun,type="marginal") {
  times=sort(unique(dat$death_time[dat$death==1]))
  BaseLambda12=basehaz(cox12,centered=F)
  BaseLambda12_fun=stepfun(BaseLambda12$time,c(0,BaseLambda12$hazard))
  BaseLambda13=basehaz(cox13,centered=F)
  BaseLambda13_fun=stepfun(BaseLambda13$time,c(0,BaseLambda13$hazard))
  BaseLambda23=basehaz(cox23,centered=F)
  BaseLambda23_fun=stepfun(BaseLambda23$time,c(0,BaseLambda23$hazard))
  baselambda12=BaseLambda12_fun(times)
  baselambda13=BaseLambda13_fun(times)
  baselambda23=BaseLambda23_fun(times)
  if(type=="marginal") {
    explp12=matrix(1,nrow=length(unique(dat$id)),ncol=length(times))
    explp13=matrix(1,nrow=length(unique(dat$id)),ncol=length(times))
    explp23=matrix(1,nrow=length(unique(dat$id)),ncol=length(times))
  } else {
    explp12=do.call(rbind,lapply(split(dat,dat$id),explp_calc,cox=cox12,times=times,ttfun=ttfun,type=type))
    explp13=do.call(rbind,lapply(split(dat,dat$id),explp_calc,cox=cox13,times=times,ttfun=ttfun,type=type))
    explp23=do.call(rbind,lapply(split(dat,dat$id),explp_calc,cox=cox23,times=times,ttfun=ttfun,type=type))
  }

  return(list(baselambda12=baselambda12,baselambda13=baselambda13,baselambda23=baselambda23,
              explp12=explp12,explp13=explp13,explp23=explp23))
}

hum_step=function(eval_times,VUS_result) {
  conc=VUS_result[[2]]
  comp=VUS_result[[3]]
  time=VUS_result[[4]]
  conc=conc[order(time)]
  comp=comp[order(time)]
  time=time[order(time)]
  vals=conc/comp
  vals_time=time
  vals_fun=stepfun(vals_time,c(0,vals))
  return(vals_fun(eval_times))
}

ivus=function(times,VUS_result,cox12_marginal,cox13_marginal,cox23_marginal,ttfun) {
  modvecs=get_modvecs(dat,cox12_marginal,cox13_marginal,cox23_marginal,ttfun=ttfun,type="marginal")
  result=P(times,modvecs$baselambda12,modvecs$explp12,
           modvecs$baselambda13,modvecs$explp13,modvecs$baselambda23,modvecs$explp23)[[1]]
  indices=which(times %in% VUS_result[[4]])
  P1=result[1,indices]
  P2=result[2,indices]
  P3=result[3,indices]
  P3_temp=result[3,indices]-c(0,result[3,indices][-length(indices)])
  
  return(sum(VUS_result[[1]]*P1*P2*P3_temp,na.rm=T)/sum(P1*P2*P3_temp,na.rm=T))
}

Nsim=1000
n=100
tmax=30
eval_times=1:tmax
c_indices_A=matrix(NA,Nsim,3)
c_indices_B=matrix(NA,Nsim,3)
HUM_ID_A=matrix(NA,Nsim,length(eval_times))
HUM_CD_A=matrix(NA,Nsim,length(eval_times))
HUM_ID_B=matrix(NA,Nsim,length(eval_times))
HUM_CD_B=matrix(NA,Nsim,length(eval_times))
IVUS_A=vector("numeric")
IVUS_B=vector("numeric")
time_vals=vector("numeric")
coef_mat=matrix(NA,Nsim,4)
for(i in 1:Nsim) {
  show(i)
  dat=get_data(n,tmax)
  dat_trans23=dat[dat$ill==1,]
  
  ttfun=function(t) {return(t > tmax/2)}
  
  dat_tt_ill=survSplit(Surv(ill_time,ill)~.,dat,cut=dat$ill_time)
  dat_tt_ill$Z2_tt=dat_tt_ill$Z2*ttfun(dat_tt_ill$ill_time)
  dat_tt_death=survSplit(Surv(ill_time,death1)~.,dat,cut=dat$ill_time)
  dat_tt_death$Z2_tt=dat_tt_death$Z2*ttfun(dat_tt_death$ill_time)
  dat_tt_trans23=survSplit(Surv(death_time,death)~.,dat_trans23,cut=dat_trans23$death_time)
  dat_tt_trans23$Z2_tt=dat_tt_trans23$Z2*ttfun(dat_tt_trans23$death_time)
  dat_tt_trans23$tstart=dat_tt_trans23$tstart+dat_tt_trans23$ill_time
  dat_tt_trans23$death_time=dat_tt_trans23$tstart+dat_tt_trans23$death_time
  
  cox12_A=coxph(Surv(tstart,ill_time,ill)~Z1,data=dat_tt_ill,timefix=F)
  cox13_A=coxph(Surv(tstart,death_time,death1)~Z1,data=dat_tt_death,timefix=F)
  cox23_A=coxph(Surv(tstart,death_time,death)~Z1,data=dat_tt_trans23,timefix=F)
  
  cox12_B=coxph(Surv(tstart,ill_time,ill)~Z1+Z2+Z2_tt,data=dat_tt_ill,timefix=F)
  cox13_B=coxph(Surv(tstart,ill_time,death1)~Z1+Z2+Z2_tt,data=dat_tt_death,timefix=F)
  cox23_B=coxph(Surv(tstart,death_time,death)~Z1+Z2+Z2_tt,data=dat_tt_trans23,timefix=F)
  
  cox12_marginal=coxph(Surv(tstart,ill_time,ill)~1,data=dat_tt_ill,timefix=F)
  cox13_marginal=coxph(Surv(tstart,ill_time,death1)~1,data=dat_tt_death,timefix=F)
  cox23_marginal=coxph(Surv(tstart,death_time,death)~1,data=dat_tt_trans23,timefix=F)
  
  c_indices_A[i,]=c(concordance(cox13_A)$concordance,concordance(cox12_A)$concordance,concordance(cox23_A)$concordance)
  c_indices_B[i,]=c(concordance(cox13_B)$concordance,concordance(cox12_B)$concordance,concordance(cox23_B)$concordance)
  
  modvecs_A=get_modvecs(dat_tt_death,cox12_A,cox13_A,cox23_A,ttfun=ttfun,type="A")
  modvecs_B=get_modvecs(dat_tt_death,cox12_B,cox13_B,cox23_B,ttfun=ttfun,type="B")
  
  times=sort(unique(dat$death_time[dat$death==1]))
  
  ill_times=dat$ill_time
  ill_times[dat$ill_time==dat$death_time & dat$ill==0]=Inf
  death_times=dat$death_time
  death_times[dat$death==0]=Inf
  censoring_times=rep(Inf,length(death_times))
  
  VUS_ID_A=VUS_ID(times,ill_times,death_times,censoring_times,modvecs_A$baselambda12,modvecs_A$explp12,
                  modvecs_A$baselambda13,modvecs_A$explp13,modvecs_A$baselambda23,modvecs_A$explp23)
  VUS_ID_B=VUS_ID(times,ill_times,death_times,censoring_times,modvecs_B$baselambda12,modvecs_B$explp12,
                  modvecs_B$baselambda13,modvecs_B$explp13,modvecs_B$baselambda23,modvecs_B$explp23)
  VUS_CD_A=VUS_CD(times,ill_times,death_times,modvecs_CD_A$baselambda12,modvecs_CD_A$explp12,
                modvecs_CD_A$baselambda13,modvecs_CD_A$explp13,modvecs_CD_A$baselambda23,modvecs_CD_A$explp23)
  VUS_CD_B=VUS_CD(times,ill_times,death_times,modvecs_CD_B$baselambda12,modvecs_CD_B$explp12,
                modvecs_CD_B$baselambda13,modvecs_CD_B$explp13,modvecs_CD_B$baselambda23,modvecs_CD_B$explp23)
  
  HUM_ID_A[i,]=hum_step(eval_times,VUS_ID_A)
  HUM_ID_B[i,]=hum_step(eval_times,VUS_ID_B)
  HUM_CD_A[i,]=hum_step(eval_times,VUS_CD_A)
  HUM_CD_B[i,]=hum_step(eval_times,VUS_CD_B)
  
  time_vals[i]=max(VUS_ID_A[[4]])
  
  IVUS_A[i]=ivus(times,VUS_ID_A,cox12_marginal,cox13_marginal,cox23_marginal,ttfun)
  IVUS_B[i]=ivus(times,VUS_ID_B,cox12_marginal,cox13_marginal,cox23_marginal,ttfun)
}
HUM_ID_A=HUM_ID_A[,-1]
HUM_ID_B=HUM_ID_B[,-1]
HUM_CD_A=HUM_CD_A[,-1]
HUM_CD_B=HUM_CD_B[,-1]
apply(c_indices_A,2,mean)
apply(c_indices_B,2,mean)
apply(HUM_ID_A,2,mean)
apply(HUM_CD_A,2,mean)
apply(HUM_ID_B,2,mean)
apply(HUM_CD_B,2,mean)
mean(IVUS_A)
mean(IVUS_B)
apply(coef_mat,2,mean)

plot(eval_times[-1],apply(HUM_ID_B,2,mean),type="l",ylim=c(0,0.75))
lines(eval_times[-1],apply(HUM_ID_A,2,mean),col=2)



