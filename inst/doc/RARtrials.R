## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(RARtrials)

## ----eval=FALSE---------------------------------------------------------------
#  #Example: RPTW(1,1) with the first 1 represents the initial number of balls for
#  #each treatment group in the urn and the second 1 represents the number of balls
#  #added to the urn when result of each participant becomes available.
#  #The function call below selects the cut-off value to control the type I error.
#  set.seed(12345)
#  sim1a<-lapply(1:5000, function(x){
#    sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
#        h=c(0.5,0.5),N2=192,side='upper')})
#  sum(sapply(sim1a, "[[", 2)>1.988,na.rm=T)/5000 #0.025
#  #Select a cut-off value to attain type I error 0.025 using the sample size 192
#  
#  sim1b<-lapply(1:5000, function(x){
#    sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
#        h=c(0.5,0.7),N2=192,side='upper',Z=1.988)})
#  sum(sapply(sim1b, "[[", 1)==1)/5000
#  #Using the selected cut-off value 1.988, we obtain power 0.7938, which is close to 0.8
#  
#  #Example: RPTW(1,1) with the first 1 represents the initial number of balls for
#  #each treatment group in the urn and the second 1 represents the number of balls
#  #added to the urn when result of each participant becomes available.
#  #Directly using asymptotic chi-square test which is equivalent to Z test statistics
#  set.seed(12345)
#  sim1<-lapply(1:5000, function(x){
#    sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
#        h=c(0.5,0.7),N2=192,side='upper')})
#  sum(sapply(sim1, "[[", 1)==1)/5000
#  #Using standard Z test statistics from normal distribution, we obtain power of 0.8038

## -----------------------------------------------------------------------------
#Example: Doubly adaptive biased coin design with five arms targeting 
#RSIHR allocation using minimal variance strategy with return of allocation 
#probabilities before applying Hu \& Zhang's formula.
dabcd_min_var(NN=c(20,23,18,25,27),Ntotal1=c(54,65,72,60,80),armn=5,type='RSIHR',
                dabcd=FALSE,gamma=2)
#The function call return values:0.3014955 0.1942672 0.0960814 0.2887962 0.1193597

#Doubly adaptive biased coin design with five arms targeting RSIHR 
#allocation using minimal variance strategy with return of allocation 
#probabilities after applying Hu \& Zhang's formula.
dabcd_min_var(NN=c(20,23,18,25,27),Ntotal1=c(54,65,72,60,80),armn=5,type='RSIHR',
                dabcd=TRUE,gamma=2)
#The function call return values:0.2076180 0.2029166 0.1717949 0.2195535 0.1981169

## -----------------------------------------------------------------------------
#Example: Doubly adaptive biased coin design with three arms targeting 
#RSIHR allocation using maximal power strategy with return of allocation
#probabilities before applying Hu \& Zhang's formula.
dabcd_max_power(NN=c(20,60,60),Ntotal1=c(86,90,90),armn=3,BB=0.1, type='Neyman',dabcd=FALSE)
#The function call return values:0.4741802 0.2629099 0.2629099
dabcd_max_power(NN=c(20,33,34),Ntotal1=c(86,78,90),armn=3,BB=0.1, type='Neyman',dabcd=FALSE)
#The function call return values:0.4433424 0.4566576 0.1000000

#Doubly adaptive biased coin design with three arms targeting RSIHR 
#allocation using maximal power strategy with return of allocation
#probabilities after applying Hu \& Zhang's formula.
dabcd_max_power(NN=c(20,60,60),Ntotal1=c(86,90,90),armn=3,BB=0.1, type='Neyman',dabcd=TRUE)
#The function call return values:0.7626214 0.1186893 0.1186893
dabcd_max_power(NN=c(20,33,34),Ntotal1=c(86,78,90),armn=3,BB=0.1, type='Neyman',dabcd=TRUE)
#The function call return values:0.427536837 0.567983270 0.004479893

## ----eval=FALSE---------------------------------------------------------------
#  #Example: Generalized RSIHR optimal allocation with unknown variances in three-armed trials
#  set.seed(12345)
#  #Under the null hypothesis
#  sim2a<-lapply(1:5000,function(x){sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,
#  TimeToOutcome=expression(rnorm( length( vStartTime ),30,3)), enrollrate=0.1,N1=12,N2=132,
#  armn=3,mean=c(9.1/100,9.1/100,9.1/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#  armlabel = c(1,2,3),cc=mean(c(9.1/100,9.1/100,9.1/100)),side='lower')})
#  h0decision<-t(sapply(sim2a, "[[", 1))
#  sum(h0decision[,1]==1|h0decision[,2]==1)/5000
#  #Attain lower one-sided type I error of 0.0218 with 5000 simulations
#  
#  #Under the alternative hypothesis
#  sim2b<-lapply(1:5000,function(x){sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,
#  TimeToOutcome=expression(rnorm( length( vStartTime ),30,3)), enrollrate=0.1,N1=12,N2=132,
#  armn=3,mean=c(9.1/100,8.47/100,8.47/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#  armlabel = c(1,2,3),cc=mean(c(9.1/100,8.47/100,8.47/100)),side='lower')})
#  h1decision<-t(sapply(sim2b, "[[", 1))
#  sum(h1decision[,1]==1)/5000
#  sum(h1decision[,2]==1)/5000
#  sum(h1decision[,1]==1|h1decision[,2]==1)/5000
#  #Marginal power of rejecting H02 is 0.8472
#  #Marginal power of rejecting H03 is 0.8432
#  #Overall power of rejecting H02 or H03 is 0.947

## -----------------------------------------------------------------------------
#### which.is.max is adapt from 'nnet' package
#### which.is.min is a rewrite from which.is.max
which.is.max <- function(x)
{
  y <- seq_along(x)[x == max(x)]
  if(length(y) > 1L) sample(y, 1L) else y
}

which.is.min <- function(x)
{
  y <- seq_along(x)[x == min(x)]
  if(length(y) > 1L) sample(y, 1L) else y
}
#Example: sample code using simulations
set.seed(12345)
#Pre-specified parameters from the Beta distribution for each arm
alpha=c(30,41,35)
beta=c(30,20,27)
#Total number of treatment groups
armn<-3
#Number of treatment groups left at current stage
armleft<-c(1,2,3) 
#Store simulation results for each arm
set.seed(12345)
result<-vector("list",length(armleft))
for (j in 1:length(armleft)) {   
  result[[j]]<- as.data.frame(rbeta(1000000,alpha[armleft[j]], 
                beta[armleft[j]]))
  colnames(result[[j]])<-sprintf("r%s",armleft[j])
}
#Expect the treatment group to have a larger treatment effect compared to the control group
#Combine results into a data frame and select the maximal value of each row
result1<-as.data.frame(do.call(cbind,result))
result1$max<-apply(result1, 1, which.is.max)
#Store results for Pr(p_k>p_{control}+\delta|data_{t-1})
theta1<-vector("list",armn)
#Store results for Pr(p_k=max{p_1,...,p_K}|data_{t-1})
pi<-vector("list",armn)
for (j in 1:length(armleft)) {
  theta1[[armleft[j]]]<-sum(result1[,j]>(result1[,1]+0.1))/1000000
  pi[[armleft[j]]]<-(sum(result1[,length(armleft)+1]==j)/1000000)
}
do.call(cbind,theta1)
#Return results: 0 0.794355 0.347715
do.call(cbind,pi)
#Return results: 0.018097 0.879338 0.102565

#Expect the treatment group to have a smaller treatment effect compared to the control group
#Combine results into a data frame and select the minimal value of each row
result1<-as.data.frame(do.call(cbind,result))
result1$max<-apply(result1, 1, which.is.min)
#Store results for Pr(p_{control}>p_k+\delta|data_{t-1})
theta1<-vector("list",armn)
#Store results for Pr(p_k=min{p_1,...,p_K}|data_{t-1})
pi<-vector("list",armn)
for (j in 1:length(armleft)) {
  theta1[[armleft[j]]]<-sum(result1[,j]<(result1[,1]-0.1))/1000000
  pi[[armleft[j]]]<-(sum(result1[,length(armleft)+1]==j)/1000000)
}
do.call(cbind,theta1)
#Return results: 0 0.001049 0.0335
do.call(cbind,pi)
#Return results: 0.755607 0.01215 0.232243

## -----------------------------------------------------------------------------

#Example: Sample code using Integration
#Expect the treatment group to have a larger treatment effect compared to the control group
#Calculate results of Pr(p_k>p_{control}+\delta|data_{t-1})
pgreater_beta(a1=alpha[1],b1=beta[1],a2=alpha[2], b2=beta[2],delta=0.1,side='upper')
pgreater_beta(a1=alpha[1],b1=beta[1],a2=alpha[3], b2=beta[3],delta=0.1,side='upper')
#Return results:  0.7951487 0.3477606

#Calculate results of Pr(p_k=max\{p_1,...,p_K\}|data_{t-1})
pmax_beta(armn=3,a2=alpha[3],b2=beta[3],a3=alpha[2],b3=beta[2],a1=alpha[1],
     b1=beta[1],side='upper')
pmax_beta(armn=3,a2=alpha[1],b2=beta[1],a3=alpha[3],b3=beta[3],a1=alpha[2],
     b1=beta[2],side='upper')
pmax_beta(armn=3,a2=alpha[1],b2=beta[1],a3=alpha[2],b3=beta[2],a1=alpha[3],
     b1=beta[3],side='upper')
#Return results: 0.01796526  0.8788907  0.1031441

#Expect the treatment group to have a smaller treatment effect compared to the control group
#Calculate results of Pr(p_{control}>p_k+\delta|data_{t-1})
pgreater_beta(a1=alpha[1],b1=beta[1],a2=alpha[2],b2=beta[2],delta=-0.1,side='lower')
pgreater_beta(a1=alpha[1],b1=beta[1],a2=alpha[3],b2=beta[3],delta=-0.1,side='lower')
#Return results:  0.001093548  0.03348547

#Calculate results of Pr(p_k=min\{p_1,...,p_K\}|data_{t-1})
pmax_beta(armn=3,a2=alpha[3],b2=beta[3],a3=alpha[2],b3=beta[2],a1=alpha[1],
     b1=beta[1],side='lower')
pmax_beta(armn=3,a2=alpha[1],b2=beta[1],a3=alpha[3],b3=beta[3],a1=alpha[2],
     b1=beta[2],side='lower')
pmax_beta(armn=3,a2=alpha[1],b2=beta[1],a3=alpha[2],b3=beta[2],a1=alpha[3],
     b1=beta[3],side='lower')
#Return results: 0.7560864  0.01230027  0.2316133

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  #Example: Select a_U by calling brar_au_binary 2000 times
#  simnull3<-lapply(1:2000,function(x){
#    set.seed(x)
#    brar_select_au_binary(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=0.9,N1=24,armn=2,
#          h=c(0.3,0.3),N2=224,tp=1,armlabel=c(1, 2),blocksize=4,alpha1=1,beta1=1,
#          alpha2=1,beta2=1,minstart=24,deltaa=-0.07,tpp=0,deltaa1=0.1,side='upper')
#  })
#  
#  #Obtain the data set of test statistics
#  simf<-list()
#  for (xx in 1:2000){
#    if (any(simnull3[[xx]][24:223,2]<0.01)){
#      simf[[xx]]<-NA
#    }  else{
#      simf[[xx]]<-simnull3[[xx]][224,2]
#    }
#  }
#  simf<-do.call(rbind,simf)
#  
#  #Ensure that around 1% of the trials stop for futility
#  sum(is.na(simf)) #20
#  #Select a_U to make sure that an overall type I error is around 0.025
#  sum(simf>0.7591,na.rm=T)/2000 #0.025
#  #The selected a_U is 0.7591.

## -----------------------------------------------------------------------------
#Example: Sample code using integration
set.seed(12345)
#Pre-specified hyper-parameters in prior distributions assumed to be the same across all 
#three groups.
para<-list(V=1/2,a=0.5,m=9.1/100,b=0.00002)
#Update hyper-parameters from the Normal-Inverse-Gamma distribution to the
#Normal-Inverse-Chi-Squared distribution.
par<-convert_gamma_to_chisq(para)
#Update hyper-parameters with some data
set.seed(123451)
y1<-rnorm(100,0.091,0.009)
par1<-update_par_nichisq(y1, par)
set.seed(123452)
y2<-rnorm(90,0.09,0.009)
par2<-update_par_nichisq(y2, par)
set.seed(123453)
y3<-rnorm(110,0.0892,0.009)
par3<-update_par_nichisq(y3, par)

#Calculate results of Pr(p_{control}>p_k+\delta|data_{t-1}) with delta=0
pgreater_NIX(par1,par2,side='lower') 
pgreater_NIX(par1,par3,side='lower') 
#Return results: 0.1959142  0.8115975
#Calculate results for Pr(p_k=min{p_1,...,p_K}|data_{t-1})
pmax_NIX(armn=3,par1=par1,par2=par2,par3=par3,side='lower') 
pmax_NIX(armn=3,par1=par2,par2=par1,par3=par3,side='lower') 
pmax_NIX(armn=3,par1=par3,par2=par2,par3=par1,side='lower') 
#Return results: 0.1801636  0.02758085  0.7922556

#Calculate results of Pr(p_k>p_{control}+\delta|data_{t-1}) with delta=0
pgreater_NIX(par1,par2,side='upper') 
pgreater_NIX(par1,par3,side='upper') 
#Return results: 0.8040858  0.1884025
#Calculate results for Pr(p_k=max\{p_1,...,p_K\}|data_{t-1})
pmax_NIX(armn=3,par1=par1,par2=par2,par3=par3,side='upper') 
pmax_NIX(armn=3,par1=par2,par2=par1,par3=par3,side='upper')
pmax_NIX(armn=3,par1=par3,par2=par2,par3=par1,side='upper') 
#Return results: 0.1876753  0.7873393  0.02498539

## ----eval=FALSE---------------------------------------------------------------
#  #Example: sample code using simulations
#  #Convert hyper-parameters to Normal-Inverse-Gamma distribution
#  convert_chisq_to_gamma<-function(cpar){
#    list(
#      m=cpar$mu,
#      V=1/cpar$kappa,
#      a=cpar$nu/2,
#      b=cpar$nu*cpar$sigsq/2
#    )
#  }
#  rnigamma<-function(nsim,par){
#    sigma2<-1/rgamma(nsim, shape=par$a,rate=par$b)
#    mu<-rnorm(nsim, par$m, sqrt(sigma2*par$V))
#    cbind(mu,sigma2)
#  }
#  set.seed(12341)
#  NIG_par1<-rnigamma(10000,convert_chisq_to_gamma(par1))
#  set.seed(12342)
#  NIG_par2<-rnigamma(10000,convert_chisq_to_gamma(par2))
#  set.seed(12343)
#  NIG_par3<-rnigamma(10000,convert_chisq_to_gamma(par3))
#  
#  #Calculate results of Pr(p_{control}>p_k+\delta|data_{t-1}) with delta=0
#  #Calculate results for Pr(p_k=min{p_1,...,p_K}|data_{t-1})
#  dat<-matrix(NA,10000,5)
#  for (i in 1:10000){
#    dat1<-rnorm(10000,NIG_par1[i,1],NIG_par1[i,2])
#    dat2<-rnorm(10000,NIG_par2[i,1],NIG_par2[i,2])
#    dat3<-rnorm(10000,NIG_par3[i,1],NIG_par3[i,2])
#    dat[i,1]<-sum(dat1>dat2)/10000
#    dat[i,2]<-sum(dat1>dat3)/10000
#    minimal<-base::pmin(dat1,dat2,dat3)
#    dat[i,3]<-sum(minimal==dat1)/10000
#    dat[i,4]<-sum(minimal==dat2)/10000
#    dat[i,5]<-sum(minimal==dat3)/10000
#  }
#  mean(dat[,1])
#  mean(dat[,2])
#  #Return results: 0.1994863 0.8113217
#  mean(dat[,3])
#  mean(dat[,4])
#  mean(dat[,5])
#  #Return results: 0.1802267 0.0285785 0.7911948
#  
#  #Calculate results of Pr(p_k>p_{control}+\delta|data_{t-1}) with delta=0
#  #Calculate results for Pr(p_k=max{p_1,...,p_K}|data_{t-1})
#  dat<-matrix(NA,10000,5)
#  for (i in 1:10000){
#    dat1<-rnorm(10000,NIG_par1[i,1],NIG_par1[i,2])
#    dat2<-rnorm(10000,NIG_par2[i,1],NIG_par2[i,2])
#    dat3<-rnorm(10000,NIG_par3[i,1],NIG_par3[i,2])
#    dat[i,1]<-sum(dat1<dat2)/10000
#    dat[i,2]<-sum(dat1<dat3)/10000
#    maximal<-base::pmax(dat1,dat2,dat3)
#    dat[i,3]<-sum(maximal==dat1)/10000
#    dat[i,4]<-sum(maximal==dat2)/10000
#    dat[i,5]<-sum(maximal==dat3)/10000
#  }
#  mean(dat[,1])
#  mean(dat[,2])
#  #Return results: 0.8005126  0.1886812
#  mean(dat[,3])
#  mean(dat[,4])
#  mean(dat[,5])
#  #Return results:  0.1910363 0.7838454 0.02511826

## ----eval=FALSE---------------------------------------------------------------
#  #Example: Select a_U by calling brar_au_binary 2000 times
#  set.seed(12345)
#  simnull4<-lapply(1:2000,function(x){
#    brar_select_au_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(
#    rnorm(length( vStartTime ),30, 3)),enrollrate=0.1, N1=192,armn=3,
#    N2=1920,tp=1,armlabel=c(1,2,3),blocksize=6,
#    mean=c((9.1/100+8.92/100+8.92/100)/3,(9.1/100+8.92/100+8.92/100)/3,
#           (9.1/100+8.92/100+8.92/100)/3),sd=c(0.009,0.009,0.009),minstart=192,
#           deltaa=c(0.00077,0.00077),tpp=1,deltaa1=c(0,0),V01=1/2,a01=0.3,m01=9/100,
#           b01=0.00001,side='lower')
#    })
#  
#  #Obtain the data set of test statistics
#  simf<-list()
#  for (xx in 1:2000){
#    if (any(simnull4[[xx]][192:1919,2]<0.01)){
#      simf[[xx]]<-NA
#    }  else{
#      simf[[xx]]<-simnull4[[xx]][1920,2]
#    }
#  }
#  simf4a<-do.call(rbind,simf)
#  
#  simf<-list()
#  for (xx in 1:2000){
#    if (any(simnull4[[xx]][192:1919,3]<0.01)){
#      simf[[xx]]<-NA
#    }  else{
#      simf[[xx]]<-simnull4[[xx]][1920,3]
#    }
#  }
#  simf4b<-do.call(rbind,simf)
#  
#  #Ensure that around 1% of the trials stop for futility
#  sum(is.na(simf4a)) #21
#  sum(is.na(simf4b)) #22
#  #Select a_U to make sure that the overall type I error is around 0.025
#  sum((simf4a[,1]>0.9836| simf4b[,1]>0.9836),na.rm=T )/2000#0.025
#  #The selected a_U is 0.9836.
#  
#  #Example to obtain the power
#  sim4<-lapply(1:1000,function(x) {
#      sim_brar_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#        length(vStartTime ),30, 3)),enrollrate=0.1,N1=192,armn=3,
#              a_U=c(0.9836,0.9836),N2=1920,tp=1,armlabel=c(1,2,3),blocksize=6,
#              mean=c(9.1/100,8.92/100,8.92/100),sd=c(0.009,0.009,0.009),
#              minstart=192,deltaa=c(0.00077,0.00077),tpp=1,deltaa1=c(0,0),
#              V01=1/2,a01=0.3,m01=9/100,b01=0.00001,side='lower')
#  
#    })
#  
#  decision<-t(sapply(sim4, "[[", 1))
#  sum(decision[,2]=='Superiorityfinal')/1000 #0.882
#  #The simulated power from 1000 simulations is 0.882.

## ----eval=FALSE---------------------------------------------------------------
#  #Example: The forward-looking Gittins Index rule applied in continuous outcomes with known
#  #variances
#  #Select cut-off values under null hypothesis for continuous outcomes with known variances.
#  #The delayed responses follow a normal distribution with a mean of 30 days, sd of 3 days
#  #and enrollment rate of 1.
#  set.seed(123452)
#  stopbound5<-lapply(1:1000,function(x){flgi_cut_off_known_var(
#  Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,
#  noRuns2=100,Tsize=120,block=8,rule='FLGI PM',prior_n=rep(1,3),prior_mean=rep(0,3),
#  mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper')})
#  stopbound5a<-do.call(rbind,stopbound5)
#  sum(stopbound5a[,1]>2.074 | stopbound5a[,2]> 2.074)/1000 #0.05
#  #The selected cut-off value is 2.074 with an overall upper one-sided type I error of 0.05.
#  #It is recommended to run more simulations to obtain a more accurate cut-off value.
#  
#  #Calculate power based on 1000 simulations
#  set.seed(123452)
#  sim5<-lapply(1:1000,function(x){sim_flgi_known_var(Gittinstype='KV',
#  df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,
#  noRuns2=100,Tsize=120,block=8,rule='FLGI PM', prior_n=rep(1,3),
#  prior_mean=rep(0,3), mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346),
#  stopbound =2.074,side='upper')})
#  h1decision<-t(sapply(sim5, "[[", 1))
#  sum(h1decision[,1]==1)/1000 #0.074
#  sum(h1decision[,2]==1)/1000 #0.244
#  sum(h1decision[,1]==1 | h1decision[,2]==1)/1000 #0.267
#  #Marginal power of rejecting H02 is 0.074
#  #Marginal power of rejecting H03 is 0.244
#  #Power of rejecting H02 or H03 is 0.267

