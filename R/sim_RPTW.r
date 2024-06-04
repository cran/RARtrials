#' @title Simulate a Trial Using Randomized Play-the-Winner Rule for Binary Endpoint
#' @description Simulate randomized play-the-winner rule in a two-armed trial with binary endpoint.
#' @details This function simulates trials using the randomized play-the-winner
#' rule under both no delay and delayed scenarios. This method is a type of urn design 
#' with the motivation to allocate more participants to the better treatment group.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_RPTW
#' @export sim_RPTW
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
#' @param na0,nb0 the initial number of balls in the urn represents the control group and the treatment group.
#' @param na1,nb1 additional number of balls represents the control group and the treatment group added to 
#' the urn after the result of each participant.
#' @param h a vector of hypothesis, for example, as c(0.1,0.1) where 0.1 stands for the success probability
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probabilities
#' for the control and a treatment group, respectively.
#' @param alphaa the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.
#' @param N2 maximal sample size for the trial.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param Z the selected cut-off value. Only specified Z when the cut-off value is selected by simulations.
#' @return \code{sim_RPTW} returns an object of class "rptw". An object of class "rptw" is a list containing 
#' final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected, 
#' Z test statistics, the simulated data set and participants accrued for each arm at the time of termination
#' of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' result with 1 stand for selected and 0 stand for not selected.
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @examples
#' #sim_RPTW with no delay responses
#' sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=0.9,na0=1,nb0=1,na1=1,nb1=1,
#' h=c(0.1,0.3),alphaa=0.025,N2=168,side='upper')
#' #sim_RPTW with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days
#' sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
#' enrollrate=0.9,na0=1,nb0=1,na1=1,nb1=1,h=c(0.1,0.3),alphaa=0.025,N2=168,side='upper')
#' @references 
#' \insertRef{Wei1978}{RARtrials}

sim_RPTW<-  function(Pats,nMax,TimeToOutcome,enrollrate,na0,nb0,na1,nb1,h,alphaa=0.025,N2,side,Z=NULL){

  start<-c(rep(1,nb0),rep(0,na0)) #1-treatment;0-control
  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome )

  data1<-matrix(NA_real_,nrow=N2,ncol=5)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  phi<-NA

  for (i in 1:N2) {
    data1[i,4]<-sample(start, 1, replace = TRUE)#1-treatment;0-control

    if (data1[i,4]==1) {
      data1[i,5]<-rbinom(1,size=1,prob=h[2]) #1 survival,0 death

      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
      if (total1>0) {

        dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[i,2])),],ncol=5)
        alive1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==1,,drop=F])#treatment alive
        alive2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==1,,drop=F])#control alive
        death1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==0,,drop=F])#treatment dead
        death2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==0,,drop=F])#control dead

      }else if (total1==0){
        death1<-0
        death2<-0
        alive1<-0
        alive2<-0
      }

      start<-c(rep(0,(alive2+death1)*na1+na0),rep(1,(death2+alive1)*nb1+nb0))

    } else if (data1[i,4]==0) {
      data1[i,5]<-rbinom(1,size=1,prob=h[1])
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
      if (total1>0) {
        dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[i,2])),],ncol=5)#data1[1:total1,,drop=F]
        alive1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==1,,drop=F])
        alive2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==1,,drop=F])
        death1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==0,,drop=F])
        death2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==0,,drop=F])

      }else if (total1==0){
        death1<-0
        death2<-0
        alive1<-0
        alive2<-0
      }

      start<-c(rep(0,(alive2+death1)*na1+na0),rep(1,(death2+alive1)*nb1+nb0))
    }
  }
  result<-data.frame(coutcome=data1[,4],doutcome=data1[,5])
  na=nrow(result[ which( result$coutcome==0),])#control
  nb=nrow(result[ which( result$coutcome==1),])
  pa=sum(nrow(result[ which( result$coutcome==0 & result$doutcome==1) , ]))/sum(nrow(result[ which( result$coutcome==0),]))
  pb=sum(nrow(result[ which( result$coutcome==1 & result$doutcome==1) , ]))/sum(nrow(result[ which( result$coutcome==1),]))
  if(is.na(pa)){pa<-0}
  if(is.na(pb)){pb<-0}
  phi<-(pb-pa)/sqrt(pa*(1-pa)/na+pb*(1-pb)/nb)

  if (side=='lower'){
    if (is.null(Z)){
      if((pnorm(phi))<=alphaa & !is.nan(phi)){
        prr1<-1 #success
      }else{
        prr1<-0
      }
    }else{
      if(phi<=Z & !is.nan(phi)){
        prr1<-1 #success
      }else{
        prr1<-0
      }
    }
    
  }else if (side=='upper'){
    if (is.null(Z)){
      if((pnorm(phi))>=(1-alphaa) & !is.nan(phi)){
        prr1<-1 #success
      }else{
        prr1<-0
      }
    }else{
      if(phi>=Z & !is.nan(phi)){
        prr1<-1 #success
      }else{
        prr1<-0
      }
    }
  }
  
  output1<-list(prr1,phi,data1,c(na,nb))
  class(output1)<-'rptw'

  return(output1)
  #return(list(prr1,phi,data1))

}


#' @export 
print.rptw<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}