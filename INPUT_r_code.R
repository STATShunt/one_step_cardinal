############################################################################################################
# READ ME
############################################################################################################

# please cite Sarkar et al. (2022) if you use substantive portions of this code for new research

# requires INPUT_data.csv and INPUT_comparison.csv in the working directory
# will produce files (in working directory) with prefix "OUTPUT_" that match (to within simulation error) the outputs from Sarkar et al. (2022)
# setting the nsim to 1000 (number of bootstrap samples used in the paper) will be slow

############################################################################################################
## data upload and set-up
############################################################################################################
cat("\014")  # clean console window
rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv) # remove all variables and start again :)

# set number of bootstrap re-samples for confidence intervals (1000 takes several minutes)
nsim<-100

# these need to be installed for the charts (use install.packages to download in not present) ...
require("tidyr")  # install.packages(tidyr)
require("ggplot2")# install.packages(ggplot2)

out_folder<-paste0(getwd(),"\\OUTPUT_")

# read in raw file
temp<-read.csv(file="INPUT_data.csv",check.names = FALSE,stringsAsFactors =TRUE , na.strings = c("N/A","#REF!","NA"," ","","-99") )  
temp_comp<-read.csv(file="INPUT_comparison.csv",check.names = FALSE,stringsAsFactors =TRUE)

str(temp)
x<-temp

###########################################################################################
##  ONE-STEP MODEL ESTIMATED
############################################################################################
# requires a data.frame called "x" that has logc, time and temp columns (numeric)
logc<-as.numeric(x$logc) ; t<-as.numeric(x$time) ; Temp<-as.numeric(x$temp)

#function (takes logc,t, Temp and vec_in as inputs) and calculates sum squared errors
foverall_optim_ROSSO<-function(vec_in){
  y0<-vec_in[1]
  y_max<-vec_in[2]+(Temp>45)*vec_in[11] #includes dummy variable for high temps
  mu_opt<-vec_in[3]
  Temp_opt<-vec_in[4]
  Temp_min<-vec_in[5]
  Temp_max<-vec_in[6]
  blam0<-vec_in[7];blam1<-vec_in[8];blam2<-vec_in[9];blam3<-vec_in[10]
  lambda<-blam0+blam1*Temp+blam2*Temp^2+blam3*Temp^3
  mu_max<- mu_opt*((Temp-Temp_max)*(Temp-Temp_min)^2 )/
    ((Temp_opt-Temp_min)*((Temp_opt-Temp_min)*(Temp-Temp_opt) -
                            (Temp_opt-Temp_max)*(Temp_opt+Temp_min-2*Temp))) # as per Rosso model
  Ft<-t+1/mu_max*log(exp(-mu_max*t)+exp(-mu_max*lambda)+exp(-mu_max*t-mu_max*lambda))
  yt<-y0+mu_max*Ft-log(1+(exp(mu_max*Ft)-1)/(exp(y_max-y0)) )
  e<-c(yt-logc)^2 # squared errors of model in-sample predictions
  sum(e) # loss function to be optimized by selection of parameter values
}

# starting values of parameters
y0<-2 ; mu_opt<-1 ;Temp_opt<-43 ; Temp_min<-5 ;Temp_max<-52
y_max<-8 ; y_max_adjust<- -5 ; blam<-c(1,0,0,0)
x0<-c(y0=y0,y_max=y_max,mu_opt=mu_opt,Temp_opt=Temp_opt,Temp_min=Temp_min,Temp_max=Temp_max
      ,blam0=blam[1],blam1=blam[2],blam2=blam[3],blam3=blam[4],y_max_adjust=y_max_adjust) 
x0_initial<-x0

# numerical optimization routine (starts at x0 and searches for parameters to minimize "sum(e)")
sim_max<-optim(x0, foverall_optim_ROSSO, method = "Nelder-Mead",   control = list(maxit = 200000,reltol=1e-16,gamma=1.5), hessian = TRUE)

#theoretical var (convert loss function to log likelihood -1/2*sum(squared error))
actual_var<-solve(sim_max$hessian/4)
sig2_hat<-(sim_max$value/2)/(nrow(x)-length(x0_initial))
estimated_sd_param<-sqrt(diag(actual_var)*sig2_hat)

out_actual<-rbind(x0,sim_max$par,estimated_sd_param)
rownames(out_actual)<-c("starting values","estimate","estimate SD")


tab_file<-paste0(out_folder,"parameter_estimates.csv")
tab_blurb<-rbind("uses Nelder-Mead with customised settings")
tab_blurb<-rbind(tab_blurb,"standard errors assume normality of residuals and are inverse of hessian (FISHER information)")
write.table(tab_blurb,tab_file,row.names=F,col.names =F ,  sep=',')
suppressWarnings( write.table(out_actual,tab_file,sep=",", row.names=T,col.names =NA, append = TRUE))

foverall_predict_ROSSO<-function(vec_in,Temp,t){
  
  y0<-vec_in[1]
  y_max<-vec_in[2]+(Temp>45)*vec_in[11]
  mu_opt<-vec_in[3]
  Temp_opt<-vec_in[4]
  Temp_min<-vec_in[5]
  Temp_max<-vec_in[6]
  blam0<-vec_in[7];blam1<-vec_in[8];blam2<-vec_in[9];blam3<-vec_in[10]
  
  lambda<-blam0+blam1*Temp+blam2*Temp^2+blam3*Temp^3

  mu_max<- ( (Temp-Temp_max)*(Temp-Temp_min)^2 )/
    (
      (Temp_opt-Temp_min)*( 
        (Temp_opt-Temp_min)*(Temp-Temp_opt) -(Temp_opt-Temp_max)*(Temp_opt+Temp_min-2*Temp)
      )
    )
  
  mu_max<-mu_max*mu_opt
  
  v<-mu_max
  
  Ft<-t+1/v*log(exp(-v*t)+exp(-mu_max*lambda)+exp(-v*t-mu_max*lambda))
  yt<-y0+mu_max*Ft-log(1+(exp(mu_max*Ft)-1)/(exp(y_max-y0)) )
  
  list(yt_hat=yt,mu_max_hat=mu_max,lambda=lambda)
  
}

foverall_predict_mu_ROSSO<-function(vec_in,Temp){
  
  y0<-vec_in[1]
  y_max<-vec_in[2]+(Temp>45)*vec_in[11]
  mu_opt<-vec_in[3]
  Temp_opt<-vec_in[4]
  Temp_min<-vec_in[5]
  Temp_max<-vec_in[6]
  blam0<-vec_in[7];blam1<-vec_in[8];blam2<-vec_in[9];blam3<-vec_in[10]
  
  lambda<-blam0+blam1*Temp+blam2*Temp^2+blam3*Temp^3

  mu_max<- ( (Temp-Temp_max)*(Temp-Temp_min)^2 )/
    (
      (Temp_opt-Temp_min)*( 
        (Temp_opt-Temp_min)*(Temp-Temp_opt) -(Temp_opt-Temp_max)*(Temp_opt+Temp_min-2*Temp)
      )
    )
  
  mu_max<-mu_max*mu_opt
  
  data.frame(mu_max_hat=mu_max,lambda=lambda,temp=Temp)
  
}

sim_max2<-optim(sim_max$par+x0*.05, foverall_optim_ROSSO, method = "BFGS",
                control = list(maxit = 1000000,reltol=1e-16), hessian = TRUE)


fitted<-sim_max$par
pred_all<-foverall_predict_ROSSO(fitted,Temp,t)

pred_rsqr<-1-var(x$logc-pred_all$yt_hat)/var(x$logc)
pred_rmse<-sqrt(mean((x$logc-pred_all$yt_hat)^2))
actual_extra<-data.frame(fit=sim_max$value, rsqr=pred_rsqr,	rmse=pred_rmse,	opt_succes=sim_max$convergence)


z<-data.frame(x,pred=pred_all$yt_hat,mu_max=pred_all$mu_max_hat,lambda=pred_all$lambda)
z$temp_factor<-paste0(z$temp," degrees")
z$res<-z$logc-z$pred


model_name<-"Baranyi and Roberts with ROSSO (one-step estimation)"

chart_pred<-ggplot()+
  geom_point(data=z,aes(x=time   ,y=logc ,colour=rep) ,alpha=.3,size=4 )+
  geom_point(data=z,aes(x=time   ,y=pred ),size=1  ,alpha=.3)+
  geom_line(data=z , aes(x=time   ,y=pred),stat="smooth", method = "loess",size=1,alpha=.3,span = .35)+
  facet_wrap(. ~temp_factor ,scales="free_x",nrow=2) +
  labs(subtitle=model_name,title="actual logc versus model predictions (smoothed line)")

chart_pred

chart_errors<-ggplot()+
  geom_point(data=z,aes(x=pred   ,y=logc ,colour=rep) ,alpha=.3,size=4 )+
  geom_line(data=data.frame(x=c(0,10),y=c(0,10)) ,aes(x=x   ,y=y) ,alpha=.3,size=1 )+
  facet_wrap(. ~temp_factor ,scales="fixed",nrow=2) +
  labs(subtitle=model_name,title="prediction versus actual logc")

temp_check<-c(4:floor(fitted[names(fitted)=="Temp_max"]),fitted[names(fitted)=="Temp_max"])
zz_dum<-foverall_predict_mu_ROSSO(fitted,temp_check)
zz<-data.frame(temp=temp_check,mu_max=zz_dum$mu_max_hat)
zz$mu_max_sqrt<-sqrt(zz$mu_max)
zz_long<-pivot_longer(zz,cols=-1)
 
chart_growth_rate<-ggplot()+ 
  geom_point(data=zz_long , aes(x=temp   ,y=value,col=name) ,alpha=.3,size=2 )+
  geom_line(data=zz_long , aes(x=temp   ,y=value,col=name),size=2,alpha=.3)+
  labs(subtitle=model_name,title="implied growth rate chart (needs error bars and comparison to literature/prior belief)")+
  xlab("temperature (degrees)")

print(chart_pred)

# paper plot
figure_1<-ggplot()+
  geom_point(data=z,aes(x=time   ,y=logc ,colour=rep) ,alpha=.3,size=4 )+
  geom_point(data=z,aes(x=time   ,y=pred ),size=1  ,alpha=.3)+
  geom_line(data=z , aes(x=time   ,y=pred),stat="smooth", method = "loess",size=1,alpha=.3,span = .35)+
  facet_wrap(. ~temp_factor ,scales="free_x",nrow=2) +
  labs(y='Observed growth (log10 CFU/g)',x='Hours')+theme(legend.position = 'bottom')+
  theme(strip.text.x = element_text(size = 12))

figure_2<-ggplot()+
  geom_point(data=z,aes(x=pred   ,y=logc ,colour=rep) ,alpha=.3,size=4 )+
  geom_line(data=data.frame(x=c(0,10),y=c(0,10)) ,aes(x=x   ,y=y) ,alpha=.3,size=1 )+
  facet_wrap(. ~temp_factor ,scales="fixed",nrow=2) +
  labs(y='Observed growth log10 CFU/g)',x='Predicted growth (log10 CFU/g)')+theme(legend.position = 'bottom')+
  theme(strip.text.x = element_text(size = 12))

csize<-c(16*2,16)*.35
cname<-"FIGURE_ONE"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(figure_1)+theme(legend.position = "none")
dev.off()

cname<-"FIGURE_TWO"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(figure_2)+theme(legend.position = "none")
dev.off()

z_long<-pivot_longer(z,-c(1:6))

res_chart_temp<-ggplot(z, aes(sample = res,alpha=.25)) + 
  stat_qq(alpha=.25,size=3) +
  facet_wrap(. ~temp_factor ,scales="free_x",nrow=2)+
  stat_qq_line(linetype = "dashed",colour="red")+
  theme(legend.position = "none")+
  labs(x="N(0,1) quantiles",y="Residual error quantiles" )+
  theme(strip.text.x = element_text(size = 12))


res_chart_all<-ggplot(z, aes(sample = res,alpha=.25)) + 
  stat_qq(alpha=.25,size=3) +
  stat_qq_line(linetype = "dashed",colour="red")+
  theme(legend.position = "none")+
  labs(x="N(0,1) quantiles",y="Residual error quantiles" )


csize<-c(16,16)*.35
cname<-"FIGURE_APPENDIX_qqnorm"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(res_chart_all)
dev.off()

csize<-c(16*2,16)*.35
cname<-"FIGURE_APPENDIX_qqnorm_by_temp"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(res_chart_temp)
dev.off()

## print out legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

fig_leg<- g_legend(figure_1+theme(legend.position = "right")) 
#csize<-c(16,16)*.5
png(filename = paste0(out_folder,"LEGEND.png"),width=csize[1],height=csize[2], units = "in",res=400)
grid::grid.newpage()
grid::grid.draw(fig_leg) 
dev.off()


current_data<-data.frame(paper="Paneer Model",zz)
comp_data<-cbind(temp_comp[,c(1,2)],mu_max=temp_comp[,3]^2,mu_max_sqrt=temp_comp[,3])
compare_all<-rbind(comp_data,current_data)
compare_all$paper<-relevel(compare_all$paper,3)

chart_compare<-ggplot()+ 
  geom_line(data=compare_all , aes(x=temp   ,y=mu_max,col=paper) ,alpha=.5,size=4 )+
  geom_point(data=compare_all[!compare_all$paper=="Paneer Model",] , aes(x=temp   ,y=mu_max,col=paper) ,alpha=1,size=3 )+
  labs(y="Maximum growth rate (log10 CFU/h)",x="Temperature (°C)") 

chart_compare_sqrt<-ggplot()+ 
  geom_line(data=compare_all , aes(x=temp   ,y=mu_max_sqrt,col=paper) ,alpha=.5,size=4 )+
  geom_point(data=compare_all[!compare_all$paper=="Paneer Model",] , aes(x=temp   ,y=mu_max_sqrt,col=paper) ,alpha=1,size=3 )+
  labs(y="Maximum growth rate sqrt(log10 CFU/h)",x="Temperature (°C)")


csize<-c(16*2,16)*.35
cname<-"FIGURE_FOUR_mu_max"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(chart_compare)
dev.off()

###########################################################################################
##  AUXILLIARY statistical summaries (for bootstrap results)
############################################################################################

fsum_grid<-function(x,do.extra=TRUE,do.autocor=FALSE) {
  outn<-list(c("table"))
  x<-as.matrix(x)
  #T<- nrow(x)
  T<-apply(!is.na(x),2,sum)
  out<-T
  
  out<-rbind( out, apply( x ,2, fsum_moment, mom=1,cent=FALSE))
  out <-rbind(out, sqrt(apply( x ,2, fsum_moment, mom=2,dfadj=1)))
  out <-rbind(out, apply( x ,2, fsum_moment, mom=3,norm=TRUE))
  out <-rbind(out, apply( x ,2, fsum_moment, mom=4,norm=TRUE))
  out<-rbind( out , T/6*out[4,]^2 + T/24*(out[5,]-3)^2 )
  pskew<-2*(1- pnorm(abs(out[4,]),0,sqrt(6/T)))
  pkurt<-2*(1-pnorm(abs(out[5,]-3),0,sqrt(24/T)))
  pjb<- 1-pchisq(out[6,],2)
  out<-rbind(out , pskew , pkurt , pjb)
  
  out.names<-c("n"  , "mean" , "sd" , "skew" , "kurt" , "JarqBera" , "pval skew" , "pval kurt" , "pval JB" )
  dimnames(out)<-list(out.names,dimnames(x)[[2]])
  
  if(do.extra){
    out.extra<- apply(x,2,function(i){c(min(i,na.rm=TRUE),max(i,na.rm=TRUE),as.vector(quantile(i,c(.025,.05,.25,.5,.75,.95,.975),na.rm=TRUE)))})
    rownames(out.extra)<-c("min", "max", "q.025","q.05", "q.25", "q.5", "q.75", "q.95", "q.975")
    out<-rbind(out,out.extra)
  }
  
  blurb<-paste0("The first three rows are sample size, mean and standard deviation."
                , "  Skewness and kurtosis are relative to a normal distribution (which has skewness 0 and kurtosis 3)."
                , "  JarqBera is a JB test statistic (a weighted average of kurtosis and skewness) used to assess normality."
                , "  The p-values for skewness, kurtosis and JB are for tests of normality (low values are evidence against normality)."
                , ifelse(do.autocor,"  The autocor(ac1) row is the first lag autocorrelation and is followed by p-values for testing the null hypothesis of zero autocorrelation (based on the Ljung-Box Q-test for one and then up to number of lags shown).","")
                , "  The minimum, maximum and quantile (q.) statistics follow in the remaining rows.")
  
  outn$table<-round(out , 9)
  outn$blurb<-blurb
  outn
}


fsum_moment<-function(x, mom , norm=FALSE , cent=TRUE , dfadj=0) { 
  x<-x[!is.na(x)]
  n<-length(x)
  if (cent ==TRUE)   { out<- sum(  (x-mean(x))^mom )/( n-dfadj) }
  else   		{ out<- sum(  (x)^mom )/( n-dfadj) }
  if (norm==TRUE)	{ sigma<- sqrt( sum( (x-mean(x))^2 )/n) ; out <- out/sigma^mom }
  out }





###########################################################################################
##  BOOTSTRAPPING
############################################################################################
foverall_predict_sim_ROSSO<-function(vec_in,Temp,t){
  
  y0<-vec_in[1]
  y_max<-vec_in[2]+(Temp>45)*vec_in[11]
  mu_opt<-vec_in[3]
  Temp_opt<-vec_in[4]
  Temp_min<-vec_in[5]
  Temp_max<-vec_in[6]
  blam0<-vec_in[7];blam1<-vec_in[8];blam2<-vec_in[9];blam3<-vec_in[10]
  
  lambda<-blam0+blam1*Temp+blam2*Temp^2+blam3*Temp^3
  
  mu_max<- ( (Temp-Temp_max)*(Temp-Temp_min)^2 )/
    (
      (Temp_opt-Temp_min)*( 
        (Temp_opt-Temp_min)*(Temp-Temp_opt) -(Temp_opt-Temp_max)*(Temp_opt+Temp_min-2*Temp)
      )
    )
  
  mu_max<-mu_max*mu_opt
  
  v<-mu_max
  
  Ft<-t+1/v*log(exp(-v*t)+exp(-mu_max*lambda)+exp(-v*t-mu_max*lambda))
  yt<-y0+mu_max*Ft-log(1+(exp(mu_max*Ft)-1)/(exp(y_max-y0)) )
  
  data.frame(yt_hat=yt,mu_max_hat=mu_max,lambda=lambda,Temp,t)
  
}

res_pool<-z$res
logc_base<-z$pred


#nsim<-100 # is set at the start of script ...
out_sim<-vector("list",nsim)
starter_values_actual<-sim_max$par

for(i in 1:nsim){
  
  #resamp residuals from model residuals
  sim_all_plus_noise<-logc_base+sample(res_pool,length(logc_base), replace=TRUE)
  t<-x$time
  Temp<-as.numeric(x$temp)
  x_boot<-data.frame(logc=sim_all_plus_noise,temp=Temp,time=t)
  
  foverall_optim_boot<-function(vec_in){
    y0<-vec_in[1]
    y_max<-vec_in[2]+(Temp>45)*vec_in[11]
    mu_opt<-vec_in[3]
    Temp_opt<-vec_in[4]
    Temp_min<-vec_in[5]
    Temp_max<-vec_in[6]
    blam0<-vec_in[7];blam1<-vec_in[8];blam2<-vec_in[9];blam3<-vec_in[10]
    
    lambda<-blam0+blam1*Temp+blam2*Temp^2+blam3*Temp^3
   
    mu_max<- ( (Temp-Temp_max)*(Temp-Temp_min)^2 )/
      (
        (Temp_opt-Temp_min)*( 
          (Temp_opt-Temp_min)*(Temp-Temp_opt) -(Temp_opt-Temp_max)*(Temp_opt+Temp_min-2*Temp)
        )
      )
    
    mu_max<-mu_max*mu_opt
    v<-mu_max
    
    Ft<-t+1/v*log(exp(-v*t)+exp(-mu_max*lambda)+exp(-v*t-mu_max*lambda))
    yt<-y0+mu_max*Ft-log(1+(exp(mu_max*Ft)-1)/(exp(y_max-y0)) )
    
    e<-c(yt-sim_all_plus_noise)^2
    sum(e)
    
  }

  # explore alt start values if simplex degenerates (try more funky adjustment)
  x0<-x0_initial

  sim_max<-optim(x0, foverall_optim_boot, method = "Nelder-Mead",
                 control = list(maxit = 1000000,reltol=1e-16,gamma=1.5)
                 ,lower = -Inf, upper = Inf, hessian = FALSE)

  out_sim[[i]]$ysim<-sim_all_plus_noise
  out_sim[[i]]$x0<-x0
  out_sim[[i]]$param<-sim_max$par
  out_sim[[i]]$opt_all<-sim_max
  out_sim[[i]]$fit<-sim_max$value
  out_sim[[i]]$x<-x_boot
  
  cat(paste("bootstrap", i,"of", nsim),"\n")
}

boot_param<-as.data.frame(t(sapply(out_sim,function(i){i$param})))

# make predictions and judge errors
pred_all<-foverall_predict_ROSSO(fitted,Temp,t)
boot_pred<-lapply(1:length(out_sim),function(i){data.frame(foverall_predict_sim_ROSSO(unlist(boot_param[i,]),as.numeric(out_sim[[i]]$x$temp),out_sim[[i]]$x$time),ysim=out_sim[[i]]$ysim)})
boot_pred<-lapply(boot_pred,function(i){data.frame(i,error=i$ysim-i$yt_hat)})
boot_ysim<-lapply(out_sim,function(i){i$ysim})

boot_success<-unlist(lapply(out_sim,function(i){i$opt_all$convergence}))==0


# model fit stuff
boot_fit<-sapply(out_sim,function(i){i$fit})
boot_rsqr<-unlist(lapply(boot_pred,function(i){1-var(i$error)/var(i$ysim)}))
boot_rmse<-unlist(lapply(boot_pred,function(i){sqrt(mean(i$error^2))}))

#plot(boot_param$Temp_max,boot_rsqr)
boot_summary<-data.frame(boot_param,fit=boot_fit,rsqr=boot_rsqr,rmse=boot_rmse,opt_success=boot_success)
pcols<-c("y0","y_max" ,"mu_opt","Temp_opt","Temp_min","Temp_max","rsqr","boot_fit","boot_rmse","opt_succes","y_max_adjust")
#boot_rsqr>.9&boot_param$Temp_max<70
out_CI<-fsum_grid(boot_summary[,colnames(boot_param)%in%pcols],do.autocor=FALSE)
boot_cor<-cor(boot_summary[,colnames(boot_param)%in%pcols],method = "pearson")
boot_cor_spearman<-cor(boot_summary[,colnames(boot_param)%in%pcols],method = "spearman")


x_sum<-out_CI
x_sum_actual<-data.frame(t(fitted[names(fitted)%in%pcols]),actual_extra)

tab_file<-paste0(out_folder,"bootstrapping.csv")
tab_blurb<-rbind(model_name,paste("parameter estimates from ",nsim,"bootstraps"),"can be used for confidence intervals etc (e.g. q.025 ot q.975)","")
write.table(tab_blurb,tab_file,row.names=F,col.names =F ,  sep=',')

suppressWarnings( write.table(rbind("actual fitted estimates",""),tab_file,sep=",", row.names=F, col.names=F, append = TRUE))
suppressWarnings( write.table(x_sum_actual,tab_file,sep=",", row.names=T,col.names =NA, append = TRUE))

suppressWarnings( write.table(rbind("","bootstraps",""),tab_file,sep=",", row.names=F, col.names=F, append = TRUE))
suppressWarnings( write.table(x_sum$table,tab_file,sep=",",col.names = NA, append = TRUE))
suppressWarnings( write.table(c("",x_sum$blur),file=tab_file , append = TRUE,  sep=',', row.names=F, col.names=F ))
suppressWarnings( write.table(rbind("","bootstrap correlations",""),tab_file,sep=",", row.names=F, col.names=F, append = TRUE))
suppressWarnings( write.table(boot_cor,tab_file,sep=",",col.names = NA, append = TRUE))


# mu analysis
boot_param_do<-boot_param[boot_success,]
temp_check_coarse<-c(4:floor(fitted[names(fitted)=="Temp_max"])-1,fitted[names(fitted)=="Temp_max"])
temp_check_end<-seq(floor(fitted[names(fitted)=="Temp_max"])-4,floor(fitted[names(fitted)=="Temp_max"]),length.out=40)
temp_check<-sort(unique(c(temp_check_coarse,temp_check_end)))

actual_mu<-data.frame(foverall_predict_mu_ROSSO(fitted,temp_check));actual_mu<-data.frame(actual_mu,mu_max_sqrt=sqrt(actual_mu$mu_max_hat))
boot_pred_mu<-lapply(1:nrow(boot_param_do),function(i){data.frame(foverall_predict_mu_ROSSO(unlist(boot_param_do[i,]),temp_check))})
boot_mu<-t(sapply(boot_pred_mu,function(i){i$mu_max_hat})) # one col for each temp in temp_check

boot_mu[boot_mu<0]<-0
boot_mu_max<-data.frame(temp=temp_check,mu_max_hat=actual_mu$mu_max_hat,t(apply(boot_mu,2,function(i){quantile(i,c(.025,.975))}) ) ,check.names = FALSE)
boot_mu_max_sqrt<-data.frame(temp=temp_check,mu_max_hat_sqrt=sqrt(actual_mu$mu_max_hat),t(apply(boot_mu,2,function(i){quantile(sqrt(i),c(.025,.975))}) ) ,check.names = FALSE)

boot_mu_max[boot_mu_max<0]<-0
boot_mu_max_sqrt[boot_mu_max_sqrt<0]<-0

# output red line
tab_file<-paste0(out_folder,"REDLINE_numbers.csv")
tab_blurb<-rbind(model_name,paste("parameter estimates are actual, 95% interval is from ",nsim,"bootstraps"),"")
write.table(tab_blurb,tab_file,row.names=F,col.names =F ,  sep=',')
suppressWarnings( write.table(boot_mu_max,tab_file,sep=",", row.names=FALSE,col.names =TRUE, append = TRUE))

boot_mu_out<-data.frame(temp=temp_check,t(boot_mu))
boot_mu_out_long<-pivot_longer(boot_mu_out,cols=-1)
boot_mu_out_long<-boot_mu_out_long[boot_mu_out_long$value<1.5&boot_mu_out_long$value>=0,]

#manual settings for publication chart (different data would need different limits ...)
mtemp<-50;mtemp_min<-5
ylim<-c(0,1.1)
xlim<-c(0,mtemp)

mtemp<-51
figure_3<-ggplot()+
  geom_line(data=boot_mu_out_long[boot_mu_out_long$temp<mtemp&boot_mu_out_long$temp>mtemp_min,], aes(x=temp,y=value,group=name),size=0.5,alpha=.05)+
  geom_line(data=actual_mu[actual_mu$temp<mtemp&actual_mu$temp>mtemp_min,], aes(x=temp,y=mu_max_hat),colour=2,size=1,alpha=1)+
  geom_line(data=boot_mu_max[boot_mu_max$temp<mtemp&boot_mu_max$temp>mtemp_min,], aes(x=temp,y=`2.5%`),colour=2,size=0.5,alpha=1,linetype = "dashed")+
  geom_line(data=boot_mu_max[boot_mu_max$temp<mtemp&boot_mu_max$temp>mtemp_min,], aes(x=temp,y=`97.5%`),colour=2,size=0.5,alpha=1,linetype = "dashed")+
  coord_cartesian(ylim = ylim,xlim=xlim)+
  scale_x_continuous(expand = expansion(add=1)) +
  scale_y_continuous(breaks = seq(0, 2, by = .2))+
  labs(y="Maximum growth rate (log10 CFU/h)",x="Temperature (°C)") 
  
csize<-c(16*2,16)*.35
cname<-"FIGURE_THREE"
png(filename = paste0(out_folder,"",cname,".png"),width=csize[1],height=csize[2], units = "in",res=400)
  print(figure_3)
dev.off()


