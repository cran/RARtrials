#' @title Simulate a Trial Using Generalized RSIHR Allocation for Continuous Endpoint with Known Variances
#' @description \code{sim_RSIHR_optimal_known_var} simulates a trial for continuous endpoints with known variances,
#' and the allocation probabilities are fixed. 
#' @details This function aims to minimize the criteria \eqn{\sum_{i=1}^{K}n_i\Psi_i}
#' with constraints \eqn{\frac{\sigma_1^2}{n_1}+\frac{\sigma_k^2}{n_k}\leq C}, where \eqn{k=2,...,K}
#' for some fixed C. It is equivalent to generalized RSIHR allocation for continuous endpoints with known variances.
#' With more than two arms the one-sided nominal level of each test is \code{alphaa} divided 
#' by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_RSIHR_optimal_known_var
#' @export sim_RSIHR_optimal_known_var
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
#' @param mean a vector of hypotheses of mean, with the first one serving as the control group.
#' @param sd a vector of hypotheses of standard deviation with the first one serving as the control group.
#' @param alphaa the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.
#' @param armlabel a vector of arm labels with an example of c(1, 2), where 1 and 2 describe
#' how each arm is labeled in a two-armed trial.
#' @param cc value in the formula of measure of treatment effectiveness, usually take the average
#' of mean responses in the hypotheses. \code{cc} is the same as C in the details.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_RSIHR_optimal_known_var} returns an object of class "RSIHRoptimal". An object of class "RSIHRoptimal" is a list containing 
#' final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected,
#' Z test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' result.
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @examples
#' #Run the function with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days under null hypothesis
#' #in a three-armed trial
#' sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(
#' rnorm(length( vStartTime ),30, 3)),enrollrate=0.9,N1=9,N2=132,armn=3,
#' mean=c(9.1/100,9.1/100,9.1/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#' armlabel=c(1,2,3),cc=mean(c(9.1/100,9.1/100,9.1/100)),side='lower')
#'
#' #Run the function with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days under alternative hypothesis
#' #in a three-armed trial
#' sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(
#' rnorm(length( vStartTime ),30, 3)),enrollrate=0.9,N1=9,N2=132,armn=3,
#' mean=c(9.1/100,8.47/100,8.47/100),sd=c(0.009,0.009,0.009),alphaa=0.025,
#' armlabel=c(1,2,3),cc=mean(c(9.1/100,8.47/100,8.47/100)),side='lower')
#' @references 
#' \insertRef{Biswas2011}{RARtrials}

sim_RSIHR_optimal_known_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,N2,armn,mean,sd,alphaa=0.025,armlabel,cc,side){

   popdat<-pop(Pats,nMax,enrollrate)
   data1<-startfun1(popdat,TimeToOutcome,blocksize=armn,N1=N1,armn=armn,armlabel=armlabel,N2=N2,mean=mean,sd=sd)


  for (i in N1:(N2-1)){

    if (i<N2){
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
    }else if (i==N2){
      total1<-N2
    }

    if (total1>0 ){

      data2<-matrix(data1[which(as.numeric(data1[1:i,3])<=as.numeric(data1[i,2])),],ncol=5)
      meann<-rep(NA_real_,armn)
      phi<-vector("list",armn)

      for (j in 1:armn){
        if (length(data2[data2[,4]==j,5])>0){

          meann[j]<-mean(c(data2[data2[,4]==j,5]))

          if (j==1){
            phi[[j]]<-((sd[j])^2)/pnorm(cc,mean=meann[j],sd=sd[j],lower.tail = F)
          } else {
            phi[[j]]<-((sd[j])^2)*pnorm(cc,mean=meann[j],sd=sd[j],lower.tail = F)
          }

        }else{
          meann[j]<-NA
        }
      }

      if (!any(is.na(meann)) && length(meann)==armn ){

        phi1<-vector("list",armn)
        for (j in 1:armn){
          phi1[[j]]<-(sd[j])^2
        }

        pho1<-unlist(phi1)
        pho11<-unlist(phi)
        pho2<-sqrt(pho11[1])*sqrt(sum(pho11[2:armn]))+sum(pho1[2:armn])
        altprob<-vector("list",armn)
        altprob[[1]]<-sqrt(pho11[1])*sqrt(sum(pho11[2:armn]))/pho2
        for (j in 2:armn){
          altprob[[j]]<-pho1[j]/pho2
        }

        altprob1<-unlist(altprob)
        assign1<-sample(1:armn,size =1, prob = altprob1,replace = TRUE)#sample(pho1,N2,replace=TRUE)

      }

    }

     if (total1==0 | (total1>0 && any(is.na(meann)))  ){
       data2<-matrix(0,nrow=1,ncol=5)
       assign1<-sample(1:armn,size =1, prob = rep(1/armn,armn),replace = TRUE)

     }

    data1[i+1,4]=assign1
    data1[i+1,5]<-rnorm(1,mean[assign1],sd[assign1])
  }

  p<-matrix(NA,armn,3)
  for (m in 1:armn) {
    p[m,1]<-mean(as.numeric(data1[data1[,4]==m,5]))
    p[m,2]<-sd[m]
    p[m,3]<-length(as.numeric(data1[data1[,4]==m,5]))
  }

  pr<-vector("list",armn-1)
  prr1<-vector("list",armn-1)


  for (l in 2:armn) {
    pr[[l-1]]<-((p[l,1]-p[1,1]))/sqrt((p[l,2]*p[l,2]/p[l,3])+(p[1,2]*p[1,2]/p[1,3]))
    if (side=='lower'){
      if((pnorm(pr[[l-1]]))<=(alphaa/(armn-1))){
        prr1[[l-1]]<-1 #success
      }else{
        prr1[[l-1]]<-0
      }
    }else if (side=='upper'){
      if((pnorm(pr[[l-1]]))>=(1-alphaa/(armn-1))){
        prr1[[l-1]]<-1 #success
      }else{
        prr1[[l-1]]<-0
      }
    }
  }


  pr10<-do.call(cbind,prr1)
  zz<-do.call(cbind,pr)

 # return(list(pr10,zz,data1))
  output1<-list(pr10,zz,data1,p[,3])
  class(output1)<-'RSIHRoptimal'
  
  return(output1)
}



#' @export  
print.RSIHRoptimal<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}