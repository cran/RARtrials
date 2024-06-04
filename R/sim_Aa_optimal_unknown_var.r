#' @title Simulate a Trial Using Aa-Optimal Allocation for Continuous Endpoint with Unknown Variances
#' @description \code{sim_Aa_optimal_unknown_var} simulates a trial for continuous endpoints with unknown variances,
#' and the allocation probabilities change based on results of accumulated participants in the trial.
#' @details This function aims to minimize the criteria \eqn{tr[A^TM^{-1}(\rho)A]}
#' and minimize the overall variance of pairwise comparisons. It is analogous to Neyman
#' allocation, favoring a higher allocation ratio to the control group. With more than two treatment 
#' groups the one-sided nominal level of each test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_Aa_optimal_unknown_var
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export sim_Aa_optimal_unknown_var
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
#' @param mean a vector of hypotheses of mean for all arms in the trial,
#' with the first one serving as the control group.
#' @param sd a vector of hypotheses of standard deviation for allarms in the trial,
#' with the first one serving as the control group.
#' @param alphaa the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.
#' @param armlabel a vector of arm labels with an example of c(1, 2), where 1 and 2 describes
#' how each arm is labeled in a two-armed trial.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_Aa_optimal_known_var} returns an object of class "Aaoptimal". An object of class "Aaoptimal" is a list containing 
#' final decision based on the T test statistics with 1 stands for selected and 0 stands for not selected,
#' T test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' result.
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom stats qt
#' @examples
#' #Run the function with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days under null hypothesis
#' #in a three-armed trial
#' sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(
#' rnorm( length( vStartTime ),30, 3)),enrollrate=0.1,N1=12,N2=132,armn=3,
#' mean=c(9.1/100,9.1/100,9.1/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#' armlabel = c(1,2,3),side='upper')
#'
#' #Run the function with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days under alternative hypothesis
#' #in a three-armed trial
#' sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(
#' rnorm( length( vStartTime ),30, 3)),enrollrate=0.1,N1=12,N2=132,armn=3,
#' mean=c(9.1/100,9.28/100,9.28/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#' armlabel = c(1,2,3),side='upper')
#' @references 
#' \insertRef{Oleksandr2013}{RARtrials}

sim_Aa_optimal_unknown_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,N2,armn,mean,sd,alphaa=0.025,armlabel,side){

  popdat<-pop(Pats,nMax,enrollrate)
  data1<-startfun1(popdat,TimeToOutcome=TimeToOutcome,blocksize=2*armn,N1=N1,armn=armn,armlabel=armlabel,N2=N2,mean=mean,sd=sd)

  for (i in N1:(N2-1)){

    if (i<N2){
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
    }else if (i==N2){
      total1<-N2
    }

    if (total1>0){
      data2<-as.data.frame(matrix(data1[which(as.numeric(data1[1:i,3])<=as.numeric(data1[i,2])),],ncol=5))
      sdd<-with(data2, aggregate(V5 ~ V4, FUN = 'sd'))[,2]
      if ( all(!is.na(sdd)) & length(sdd)==armn){
        pho<-vector("list",armn)
        pho[[1]]<-sdd[1]*sqrt(armn-1)/(sdd[1]*sqrt(armn-1)+sum(sdd[2:armn]))
        for (j in 2:armn){
          pho[[j]]<-sdd[j]/(sdd[1]*sqrt(armn-1)+sum(sdd[2:armn]))
        }

        pho1<-unlist(pho)

        assign1<-sample(1:armn,size =1, prob = pho1,replace = TRUE)#sample(pho1,N2,replace=TRUE)
      } else {
        assign1<-sample(1:armn,size =1, prob = rep(1/armn,armn),replace = TRUE)
      }
    }

    if (total1==0 ){
      assign1<-sample(1:armn,size =1, prob = rep(1/armn,armn),replace = TRUE)
    }

    data1[i+1,4]=assign1
    data1[i+1,5]<-rnorm(1,mean[assign1],sd[assign1])
  }


  p<-matrix(NA,armn,3)
  for (m in 1:armn) {
    p[m,1]<-mean(as.numeric(data1[data1[,4]==m,5]))
    p[m,2]<-sd(as.numeric(data1[data1[,4]==m,5]))
    p[m,3]<-length(as.numeric(data1[data1[,4]==m,5]))
  }


  pr<-vector("list",armn-1)
  zz<-vector("list",armn-1)
  for (l in 2:armn) {
    zz[[l-1]]<-((p[l,1]-p[1,1]))/sqrt((p[l,2]*p[l,2]/p[l,3])+(p[1,2]*p[1,2]/p[1,3]))
    if (side=='lower'){
      if (zz[[l-1]]<=qt(alphaa/(armn-1),df=p[l,3]+p[1,3]-2)){
        pr[[l-1]]<-1 #success
      }else{
        pr[[l-1]]<-0
      }
    }else if (side=='upper'){
      if (zz[[l-1]]>=qt(1-alphaa/(armn-1),df=p[l,3]+p[1,3]-2)){
        pr[[l-1]]<-1 #success
      }else{
        pr[[l-1]]<-0
      }
    }
  }
  pr1<-do.call(cbind,pr)
  zz1<-do.call(cbind,zz)
 # return(list(pr1,zz1,data1))
  output1<-list(pr1,zz1,data1,p[,3])
  class(output1)<-'Aaoptimal'
  
  return(output1)
}

#' @export 
print.Aaoptimal<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}