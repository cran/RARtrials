#' @title Simulate a Trial Using Doubly Adaptive Biased Coin Design with Minmial Variance Strategy for Binary Endpoint
#' @description \code{sim_dabcd_min_var} can be used for doubly adaptive biased coin design with minimal variance
#' strategy for binary outcomes, targeting generalized Neyman allocation and generalized RSIHR allocation
#' with 2-5 arms. 
#' @details The output of this function is based on Hu \code{\&} Zhang's formula \insertCite{Hu2004}{RARtrials}.
#' With more than two arms the one-sided nominal level of each
#' test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_dabcd_min_var
#' @export sim_dabcd_min_var
#' @param Pats the number of patients accrued within a certain time frame indicates the
#' count of individuals who have been affected by the disease during that specific period,
#' for example, a month or a day. If this number is 10, it represents that
#' 10 people have got the disease within the specified time frame.
#' @param nMax the assumed maximum accrued number of patients with the disease in the population, this number
#' should be chosen carefully to ensure a sufficient number of patients are simulated,
#' especially when considering the delay mechanism.
#' @param TimeToOutcome the distribution of delayed response times or a fixed delay time for responses.
#' The delayed time could be a month, a week or any other time frame. When the unit changes,
#' the number of TimeToOutcome should also change. It can be in the format
#' of expression(rnorm( length( vStartTime ),30, 3)), representing delayed responses
#' with a normal distribution, where the mean is 30 days and the standard deviation is 3 days.
#' @param enrollrate probability that patients in the population can enroll in the trial.
#' This parameter is related to the number of people who have been affected by the disease in the population,
#' following an exponential distribution.
#' @param N1 number of participants with equal randomization in the 'initialization' period.
#' Recommend using 10 percent of the total sample size.
#' @param N2 maximal sample size for the trial.
#' @param armn number of total arms in the trial.
#' @param armlabel a vector of arm labels with an example of c(1, 2), where 1 and 2 describe
#' how each arm is labeled in a two-armed trial.
#' @param h a vector of success probabilities in hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the success probability
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probabilities
#' for the control and the treatment group, respectively.
#' @param type allocation type, with choices from 'RSIHR' and 'Neyman'.
#' @param gamma tuning parameter in Hu & Zhang's formula (\insertCite{Hu2004}{RARtrials}). When dabcd=0, this parameter does not need
#' to be specified. Default value is set to 2.
#' @param alphaa the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_dabcd_min_var} returns an object of class "dabcd". An object of class "dabcd" is a list containing 
#' final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected,
#' Z test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' result.
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @examples
#' sim_dabcd_min_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length( vStartTime ),30, 3)),enrollrate=0.9,N1=30,N2=300,armn=3,
#' armlabel=c(1,2,3),h=c(0.2,0.3,0.2),type='Neyman',
#' side='upper')
#' sim_dabcd_min_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length( vStartTime ),60, 3)),enrollrate=0.1,N1=50,N2=500,armn=3,
#' armlabel=c(1,2,3),h=c(0.2,0.3,0.3),type='RSIHR',
#' side='lower')
#' @references 
#' \insertRef{Hu2004}{RARtrials}

sim_dabcd_min_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,N2,armn,armlabel,h,type,gamma=2,alphaa=0.025,side){

  alr<-rep(NA_real_,armn)
  ncp<-rep(0,1)
  pr1<-rep(0,armn)

  popdat<-pop(Pats,nMax,enrollrate)
  data1<-startfun(popdat,TimeToOutcome=TimeToOutcome,blocksize=2*armn,N1=N1,armn=armn,armlabel=armlabel,N2=N2,h=h)
  low <- rep(0,armn)
  high <-rep(1,armn)

  for (i in N1:N2) {
    if (i>N1){
      assigna<-sample(1:armn,size = 1,prob = alr)
      data1[i,4]<- assigna
      data1[i,5]<-rbinom(1,size=1,prob=h[assigna])
    }
    if (i<N2){
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
    }else if (i==N2){
      total1<-N2
    }

    if (total1>0){
      data2<-matrix(data1[which(as.numeric(data1[1:i,3])<=as.numeric(data1[i,2])),],ncol=5)
    }else if (total1==0){
      data2<-matrix(0,nrow=1,ncol=5)
    }

    dummy1<-matrix(NA,armn,2)
    for (m in 1:armn) {
      dummy1[m,1]<-sum(as.numeric(data2[data2[,4] %in% m,5]))
      dummy1[m,2]<-nrow(data2[data2[,4]==m,,drop=F])
    }
    NN<-dummy1[,1]+1
    Ntotal1<-dummy1[,2]+2
    p<-cbind(p=unname(unlist(NN/Ntotal1)),arm=1:armn)
    rho<-rep(NA,armn)
    rho1<-rep(NA,armn)

  if ( type=='Neyman'){
    for (k in 1:armn ){
      rho1[k]<-sqrt(p[k,1]*(1-p[k,1]))
    }
    for (k in 1:armn ){
      rho[k]<-rho1[k]/(sum(rho1))
    }
  }
  if (type=='RSIHR'){
    for (k in 1:armn ){
      rho1[k]<-sqrt(p[k,1])
    }
    for (k in 1:armn ){
      rho[k]<-rho1[k]/(sum(rho1))
    }
  }

    alr<-rep(NA,armn)
    phi<-rep(NA,armn)

    if ( any(dummy1[,2]==0) ){
      phi<-rep(1/armn,armn)
    }else{
      for (k in 1:armn){
        phi[k]<-rho[k]*((rho[k]/(dummy1[k,2]/(sum(dummy1[,2]))))^gamma)
      }
    }
    for (kk in 1:armn){
      alr[kk]<-phi[kk]/sum(phi)
    }
  }


  pr<-vector("list",armn-1)
  phi<-vector("list",armn-1)
  dummy1<-matrix(NA,armn,2)
  for (m in 1:armn) {
    dummy1[m,1]<-sum(as.numeric(data1[data1[,4] %in% m,5]))
    dummy1[m,2]<-nrow(data1[data1[,4]==m,])
  }

  NN<-dummy1[,1]
  Ntotal1<-dummy1[,2]
  p<-cbind(p=unname(unlist(NN/Ntotal1)),arm=1:armn)

  for (l in 2:armn) {

    phi[[l-1]]<-(p[l,1]-p[1,1])/sqrt((p[l,1]*(1-p[l,1])/Ntotal1[l])+(p[1,1]*(1-p[1,1])/Ntotal1[1]))

    if (side=='upper'){
       if (phi[[l-1]]>=qnorm(1-alphaa/(armn-1))){
         pr[[l-1]]<-1 #success
       }else{
         pr[[l-1]]<-0
       }
    }else if (side=='lower'){
      if (phi[[l-1]]<=qnorm(alphaa/(armn-1))){
        pr[[l-1]]<-1 #success
      }else{
        pr[[l-1]]<-0
      }
    }
  }
  pr1<-do.call(cbind,pr)
  phi1<-do.call(cbind,phi)
 # return(list(pr1,phi1,data1))
  output1<-list(pr1,phi1,data1,Ntotal1)
  class(output1)<-'dabcd'
  
  return(output1)
}

#' @export 
print.dabcd<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}