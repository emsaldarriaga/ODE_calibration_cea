#Title: Simple calibration of a compartmental model and applications in economic analysis
#Author: Enrique M. Saldarriaga
#Meeting: Prevention Effectiveness Fellowship - Workshop on Simulation Methods
#Date: Oct 5, 2023

#Aim: To calibrate an SEIR model using direct search and present how the results can be leveraged in economic analysis 

#Require
pkgs = c("deSolve","dplyr","ggplot2","poissondisc","lhs","plotly","invgamma")
sapply(c(pkgs),require, character.only=T)
rm(pkgs)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Calibration of an SEIR model####
##Create model function####
SEVIR.model <- function (t, x, params) {
  
  week = x[1]
  S = x[2] #susceptible
  E = x[3] #exposed, pre-infectious
  I = x[4] #symptomatic, infectious
  R = x[5] #recovered and resistant,
  V = x[6] #vaccinated
  
  #Parameters
  N = params[1] #Initial population size; equivalent to susceptible population
  beta = params[2] #transmission parameter
  lenE = params[3] #inverse length of exposed, pre-infectious period
  lenInt = params[4] #inverse length of symptomatic, infectious period
  veff = params[5] #effectiveness of vaccine
  vr = params[6] #vaccination rate, doses per week
  
  #Probability rate of susceptibles getting infected (force of infection)
  #Denominator is mixing pool
  FOI = beta * I/N
  
  #ODE system
    d_week_dt = 1 #Weekly mode, increases by 1 week
    d_S_dt = -S * FOI - vr
    d_E_dt = (S + V*(1 - veff)) * FOI - E * lenE
    d_I_dt = E * lenE - I * lenInt
    d_R_dt = I * lenInt
    d_V_dt = vr - V*(1 - veff)*(FOI) 
    d_vr = vr
  
  return(list(c(d_week_dt,d_S_dt,d_E_dt,d_I_dt,d_R_dt,d_V_dt,d_vr)))
}

times = seq(1, 50, by=1) #Solve for 50 weeks

##Initial values of state variables
week0 = 0
S0 = 10000
E0 = 5
I0 = 0
R0 = 0
V0 = 0
vr0=0

init_values = c(week0,S0,E0,I0,R0,V0,vr0)

#Parameter values
N = 10005 #Population size (S+I)
beta = 0.7 #Transmission parameter
lenE = 7/3 #inverse length of exposed, pre-infectious period lasts 3 days
lenInt = 1/5 #inverse length of symptomatic, infectious period lasts for 5 weeks
veff = 0.45 #effectiveness of vaccine
vr=5 #vaccination rate, doses per week

parameters = c(N,beta,lenE,lenInt,veff,vr)

get_ode = function(N = 10005, beta = 0.7,lenE = 7/3, lenInt = 1/5,veff = 0.45, vr=5){
  param = c(N,beta,lenE,lenInt,veff,vr)
  r = lsoda(init_values, times, SEVIR.model, param) %>% data.frame()%>%
    setNames(c("timestep",'week',"S","E","I","R","V","CumV")) %>% mutate(N = rowSums(.[,3:7]))
  return(r)
}

results = get_ode()

ggplot(results)+geom_line(aes(x=times,y=S,col="S"),linewidth=0.8)+
  geom_line(aes(x=times,y=E,col="E"),linewidth=0.8)+
  geom_line(aes(x=times,y=I,col="I"),linewidth=0.8)+ 
  geom_line(aes(x=times,y=R,col="R"),linewidth=0.8)+
  #geom_line(aes(x=times,y=V,col="V"),linewidth=0.8)+
  labs(title="")+theme_bw()

##Simulating real-word data####
set.seed(1234)
real_I = sapply(results$I,FUN=function(x)rpois(1,x*rnorm(1,1,0.3)))

ggplot()+geom_point(aes(x=times,y=real_I,col="Real I"),size=1.8)+
  geom_path(aes(x=times,y=real_I,col="Real I"),size=0.8)+
  geom_line(aes(x=times,y=results$I,col="Sim I"),linewidth=0.8)+labs(x="Week",y="Incidence",col="",)

##Calibration: Direct/manual search####
##Target Data: I

##Free parameters: beta, lenE, lenInt; what do we know? w
#beta: Norm(0.7, sd=0.09)
#lenE: 7/gamma(shape=3,rate=1)
#lenInt: 1/gamma(shape=5,rate=2)

set.seed(1234); k = lhs::randomLHS(n = 1e3,k = 3)
set.seed(1234); p_sampling = data.frame(Index=1:1e3,beta = qnorm(k[,1],0.7,0.09),
                        lenE = 7/qgamma(k[,2],shape=3,rate=1),
                        lenInt = 1/qgamma(k[,3],shape=5,rate=1.5)) %>% 
  mutate(beta=ifelse(beta>=1,1-1e-5,beta))

#lenE: Length of pre-infectious ranges from 0.1 to 13.4 days; 
## (summary(p_sampling$lenE)/7)^-1; hist((p_sampling$lenE/7)^-1) 
#lenINt: Length of infection ranges from 0.24 to 10 weeks; 
## (summary(p_sampling$lenInt)^-1); hist(p_sampling$lenInt^-1)

##Cost function
mean(abs(results$I-real_I)) #MAE
mean((results$I-real_I)^2) #MSE
-sum(dpois(real_I,results$I,log = TRUE)) #Likelihood using Poisson distribution

##Searching algorithm
# Manual:
calibration_direct = list()
calibration_direct_s = data.frame()
tictoc::tic()
for(i in 1:nrow(k)){
  calibration_direct_s[i,"Index"] = i
  
  r_i = get_ode(beta = p_sampling[i,"beta"],
                          lenE = p_sampling[i,"lenE"],
                          lenInt = p_sampling[i,"lenInt"]) %>% 
    mutate(L_MAE = abs(I-real_I),L_MSE=(I-real_I)^2,
           L_LP = (dpois(real_I,I,log=TRUE)),
           L_LN = (dnorm(real_I,mean = I, sd = 1, log = TRUE)))
  calibration_direct[[i]] = r_i
  calibration_direct_s[i,"L_MAE"] = mean(r_i$L_MAE)
  calibration_direct_s[i,"L_MSE"] = mean(r_i$L_MSE)
  calibration_direct_s[i,"L_LP"] = -sum(r_i$L_LP)
  calibration_direct_s[i,"L_LN"] = -sum(r_i$L_LN)
}
tictoc::toc()

(ggplot(calibration_direct_s%>%arrange(L_MAE)%>%mutate(Ind=1:1e3),
       aes(x=Ind,label=Index))+
  geom_point(aes(y=L_MAE,color="L_MAE")))%>%ggplotly()

(ggplot(calibration_direct_s%>%arrange(L_MSE)%>%mutate(Ind=1:1e3),
       aes(x=Ind,label=Index))+
  geom_point(aes(y=L_MSE,color="L_MSE")))%>%ggplotly()

(ggplot(calibration_direct_s%>%arrange(L_LP)%>%mutate(Ind=1:1e3),
       aes(x=Ind,label=Index))+
  geom_point(aes(y=L_LP,color="L_LP")))%>%ggplotly()

(ggplot(calibration_direct_s%>%arrange(L_LN)%>%mutate(Ind=1:1e3),
        aes(x=Ind,label=Index))+
    geom_point(aes(y=L_LN,color="L_LN")))%>%ggplotly()

ggplot()+geom_point(aes(x=times,y=real_I,col="Real I"),size=1.8)+
  geom_line(aes(x=times,y=calibration_direct[[66]]$I,col="Sim I"),linewidth=0.8)+labs(x="Week",y="Incidence",col="")

##Posterior ranges
get_post = function(samplingdf = p_sampling,index){
  d = samplingdf %>% filter(Index%in%index)
  labss = paste0(names(d[,-1]),
                 apply(d[,-1],2,FUN=function(x)
                   paste0(": ",round(mean(x),3)," (",
                          round(quantile(x,0.025),3),", ",
                          round(quantile(x,0.975),3),")"))
                 )
  
  pbeta = ggplot(d)+geom_histogram(aes(x=beta))+
    labs(title =labss[1])
  plenE = ggplot(d)+geom_histogram(aes(x=lenE))+
    labs(title=labss[2])
  plenInt = ggplot(d)+geom_histogram(aes(x=lenInt))+
    labs(title=labss[3])
  
  print(ggpubr::ggarrange(
    pbeta,plenE,plenInt,nrow=1,ncol=3,legend = "right",common.legend = TRUE)
    )
  print(labss)
}

#The more retained realizations, the greater the chance to obtain wider confidence intervals 
ind_MAE = (sort(calibration_direct_s$L_MAE,index.return = TRUE)$ix)[1:250]
ind_MSE = (sort(calibration_direct_s$L_MSE,index.return = TRUE)$ix)[1:250]
ind_LP = (sort(calibration_direct_s$L_LP,index.return = TRUE)$ix)[1:250]
ind_LN = (sort(calibration_direct_s$L_LN,index.return = TRUE)$ix)[1:250]

get_post(index=ind_LN)

ps = ggplot()+geom_point(aes(x=times,y=real_I,col="Real I"),size=1.8)+
  geom_line()+labs(x="Week",y="Incidence",col="")

for(j in ind_MAE){
  ps = ps + geom_line(data=calibration_direct[[j]],aes(x=times,y=I,col="MAE"),
                      alpha=0.2)
}

for(j in ind_MSE){
  ps = ps + geom_line(data=calibration_direct[[j]],aes(x=times,y=I,col="MSE"),
                      alpha=0.2)
}

for(j in ind_LP){
  ps = ps + geom_line(data=calibration_direct[[j]],aes(x=times,y=I,col="LP"),
                      alpha=0.2)
}
for(j in ind_LN){
  ps = ps + geom_line(data=calibration_direct[[j]],aes(x=times,y=I,col="LN"),
                      alpha=0.2)
}
ps = ps+scale_color_manual(values=c("#0DBC98","#3636C5","#D2C72B","#053D0B","#951613"))+theme_bw(base_size=30)
# ps = ps+scale_color_manual(values=c("#ff7f00","#984ea3","#377eb8","#4daf4a","#e41a1c"))

png("Calibration.png",height = 1800,width=4/3*1800,res = 120)
ps
dev.off()

##Obtaining posterior bounds
get_post(index=Reduce(intersect,list(ind_LN,ind_MAE,ind_MSE)))

##Obtaining posterior distributions
#"beta: 0.713 (0.603, 0.843)"  "lenE: 3.55 (1.409, 10.033)"  "lenInt: 0.23 (0.158, 0.308)"

get_prms_normal = function(targetq,inits){ 
  #Target is a vector of 2 that contains the 95%CI of the desired distribution
  #Inits is a vector of the same lenght as target that contains the initial searching values; must be sensible
  get_param = function(q,target){
    d = qnorm(c(0.025,0.95),mean=q[1],sd=q[2])
    sum((d-target)^2)
  }
  op = optim(par=inits,fn=get_param,target=targetq)
  return(op$par)
}

get_prms_invgamma = function(targetq,inits){ 
  #Target is a vector of 2 that contains the 95%CI of the desired distribution
  #Inits is a vector of the same lenght as target that contains the initial searching values; must be sensible
  get_param = function(q,target){
    d = qinvgamma(c(0.025,0.95),shape=q[1],rate=q[2])
    sum((d-target)^2)
  }
  op = optim(par=inits,fn=get_param,target=targetq)
  return(op$par)
}

#beta: rnorm(1e3,0.73349,0.06657); 
get_prms_normal(c(0.603, 0.843),c(0.7, 0.1))

#lenE: rinvgamma(1e3,3.62076,11.54346); 
get_prms_invgamma(c(1.409,10.033),c(3,10))

#lenInt: rinvgamma(1e3,29.0557789,6.4044434); 
get_prms_invgamma(c(0.158,0.308),c(4,0.4))

#Drawing from posterior distributions to create realizations of model internalizing parameter uncertainty
set.seed(1234); p_posterior = data.frame(Index=1:1e3,beta = qnorm(k[,1],0.73349,0.06657),
                        lenE = qinvgamma(k[,2],3.62076,11.54346),
                        lenInt = qinvgamma(k[,3],29.0557789,6.4044434)) %>% 
  mutate(beta=ifelse(beta>=1,1-1e-5,beta))

post_realizations=vector("list",length = 6)
names(post_realizations) = c("S","E","I","R","V","CumV")
for(i in 1:nrow(k)){
  r_i = get_ode(beta = p_posterior[i,"beta"],
                lenE = p_posterior[i,"lenE"],
                lenInt = p_posterior[i,"lenInt"])
  post_realizations[["S"]] = cbind(post_realizations[["S"]],r_i$S)
  post_realizations[["E"]] = cbind(post_realizations[["E"]],r_i$E)
  post_realizations[["I"]] = cbind(post_realizations[["I"]],r_i$I)
  post_realizations[["R"]] = cbind(post_realizations[["R"]],r_i$R)
  post_realizations[["V"]] = cbind(post_realizations[["V"]],r_i$V)
  post_realizations[["CumV"]] = cbind(post_realizations[["CumV"]],r_i$CumV)
}

for(j in 1:6){
  post_realizations[[j]] = post_realizations[[j]] %>% data.frame() %>% setNames(paste0(names(post_realizations)[j],1:1000)) %>%
    cbind(Mean=apply(post_realizations[[j]],1,FUN=function(x)mean(x)),
           Median=apply(post_realizations[[j]],1,FUN=function(x)quantile(x,0.5)),
           LB=apply(post_realizations[[j]],1,FUN=function(x)quantile(x,0.025)),
           UB=apply(post_realizations[[j]],1,FUN=function(x)quantile(x,0.975)))
}

results_post_calibration = data.frame(times=times)
for(j in 1:6){
  results_post_calibration = cbind(results_post_calibration,post_realizations[[j]]%>%dplyr::select(Mean,Median,LB,UB)%>%
                                     setNames(paste0(names(post_realizations)[j],"_",colnames(.))))
}

ggplot(results_post_calibration,aes(x=times))+
  geom_line(aes(y=I_Mean,col="I"),linewidth=0.8)+
  geom_ribbon(aes(y=I_Mean,fill="I",ymin=I_LB,ymax=I_UB),alpha=0.2)+
  geom_point(aes(y=real_I,col="Real I"))+geom_path(aes(y=real_I,col="Real I"))+
  labs(title="",col="")+theme_bw()+scale_color_manual(values=c("darkblue","red"))+scale_fill_manual(values=c("darkblue"))+
  guides(fill="none")

ggplot(results_post_calibration,aes(x=times))+
  geom_line(aes(y=S_Mean,col="S"),linewidth=0.8)+
  geom_ribbon(aes(y=S_Mean,fill="S",ymin=S_LB,ymax=S_UB),alpha=0.3)+
  geom_line(aes(y=E_Mean,col="E"),linewidth=0.8)+
  geom_ribbon(aes(y=E_Mean,fill="E",ymin=E_LB,ymax=E_UB),alpha=0.3)+
  geom_line(aes(y=I_Mean,col="I"),linewidth=0.8)+
  geom_ribbon(aes(y=I_Mean,fill="I",ymin=I_LB,ymax=I_UB),alpha=0.3)+
  geom_line(aes(y=R_Mean,col="R"),linewidth=0.8)+
  geom_ribbon(aes(y=R_Mean,fill="R",ymin=R_LB,ymax=R_UB),alpha=0.3)+
  labs(title="",col="")+theme_bw()+guides(fill="none")

#Best guest draws
#"beta:  (0.603, 0.843)"  "lenE:  (1.409, 10.033)"  "lenInt:  (0.158, 0.308)"
results_best_estimation = get_ode(N = 10005, beta = 0.713, #before beta=0.7
                                  lenE = 3.55, #before 3 days (7/3), now 1.97 days 
                                  lenInt = 0.23, #before 5 weeks (1/5), now  4.347 weeks
                                  veff = 0.45, vr=5)

ggplot(results_best_estimation,aes(x=times))+
  geom_line(aes(y=I,col="I"),linewidth=0.8)+
  geom_point(aes(y=real_I,col="Real I"))+
  geom_path(aes(y=real_I,col="Real I"))+
  scale_color_manual(values=c("darkblue","red"))+
  labs(title="",col="")+theme_bw()

#Application to Economic Analysis####

## Economic analysis of healthcare innovations: improved vaccination outreach
#' Time horizon: 50 weeks (purely retrospective)
#' Discount rate: NA
#' Perspective of perceived costs and 'benefits': Built-in in cost and outcomes; healthcare sector
#' Costs -- perspective and time unit: S=0, E=10, I=20, R=0, V=8
#' Health outcomes -- perspective and time unit: S=0.92, E=0.8, I=0.76, R=0.92, V=0.93

get_cea_params = function(State=c("S","E","I","R","V"),
                      Costs = c(0,10,20,0,8), 
                      Qalys = c(0.92,0.8,0.76,0.92,0.93)){
  return(cbind(State,Costs,Qalys)%>%data.frame()%>%mutate(Costs=as.numeric(Costs),
                                                          Qalys=as.numeric(Qalys)))
}


ce_results = function(ode_results,cea_parameters){
  #WARNING: non robust to changes in names or order
  return(
    cbind(as.matrix(ode_results[,3:7])%*%(as.matrix(cea_parameters$Costs)),
        as.matrix(ode_results[,3:7])%*%(as.matrix(cea_parameters$Qalys))) %>% apply(.,2,sum)
  )
}

ce_results(results_best_estimation,get_cea_params())

results_8vr = get_ode(N = 10005, beta = 0.713,
                                  lenE = 3.55,
                                  lenInt = 0.23,
                                  veff = 0.45, vr=8)

ce_results(results_8vr, 
           get_cea_params(Costs = c(0,10,20,0,10)))

ce_results(get_ode(N = 10005, beta = 0.713, lenE = 3.55, lenInt = 0.23, veff = 0.45, vr=8), 
          get_cea_params(Costs = c(0,10,20,0,10)))

get_cea = function(x,labs=NULL,status.quo.comp=0){
  #X is a matrix of CE results, with as many rows as scenarios/interventions were analyzed
  #status.quo.comp=1 if want to include comparisson with the least effective strategy (usually this is the status quo), and 0 if not. Default is 0
  r = x %>% data.frame() %>% `rownames<-`(labs)%>%
    setNames(c("Costs","QALYs"))%>% arrange(QALYs,Costs) %>%
    mutate(dCosts = c(NA,diff(Costs)),dQALYs=c(NA,diff(QALYs)))%>%
    mutate(ICER=dCosts/dQALYs)
  rs = cbind(r,data.frame(dLE_Costs=c(NA,r[-1,1]-r[1,1]),dLE_QALYs=c(NA,r[-1,2]-r[1,2]))%>%
    mutate(LE_ICER = dLE_Costs/dLE_QALYs))
  return(if(status.quo.comp==0){r}else{rs})
}

get_cea(labs=c("Status Quo","Vr=8-higher cost"),rbind(c(872151.9, 453401.3),c( 897902.4,453437.1)))

get_cea(labs=c("Status Quo","Vr=8-higher cost","Vr=8-same cost","Vr=10-higher cost","Vr=10-same cost"),
        rbind(
          ce_results(results_best_estimation,get_cea_params()),
          ce_results(get_ode(N = 10005, beta = 0.713, lenE = 3.55, lenInt = 0.23, veff = 0.45, vr=8), 
                     get_cea_params(Costs = c(0,10,20,0,10))),
          ce_results(get_ode(N = 10005, beta = 0.713, lenE = 3.55, lenInt = 0.23, veff = 0.45, vr=8), 
                     get_cea_params()),
          ce_results(get_ode(N = 10005, beta = 0.713, lenE = 3.55, lenInt = 0.23, veff = 0.45, vr=10), 
                     get_cea_params(Costs = c(0,10,20,0,10))),
          ce_results(get_ode(N = 10005, beta = 0.713, lenE = 3.55, lenInt = 0.23, veff = 0.45, vr=10), 
                     get_cea_params())
        ),status.quo.comp = 1
        )






