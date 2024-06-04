#' @title Simulate a Trial Using Forward-Looking Gittins Index for Continuous Endpoint with Unknown Variances
#' @description Function for simulating a trial using the forward-looking Gittins Index rule and the controlled forward-looking
#' Gittins Index rule for continuous outcomes with unknown variances in trials with 2-5 arms. The prior distributions
#' follow Normal-Inverse-Gamma (NIG) (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)}) 
#' distributions and should be the same for each arm.
#' @details This function simulates a trial using the forward-looking Gittins Index rule or the
#' controlled forward-looking Gittins Index rule under both no delay and delayed scenarios.
#' The cut-off value used for \code{stopbound} is obtained by simulations using \code{flgi_stop_bound_flgi_unk_var}.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_flgi_unknown_var
#' @export sim_flgi_unknown_var
#' @param Gittinstype type of Gittins indices, should be set to 'UNKV' in this function
#' @param df discount factor which is the multiplier for loss at each additional patient in the future.
#' Available values are 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99 and 0.995. The maximal sample size can be up to 10000.
#' @param gittins user specified Gittins indices for calculation in this function. If \code{gittins} is provided,
#' \code{Gittinstype} and \code{df} should be NULL.
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
#' @param K number of total arms in the trial.
#' @param noRuns2 number of simulations for simulated allocation probabilities within each block. Default value is
#' set to 100 times, which is recommended in Villar et al., 2015.
#' @param Tsize maximal sample size for the trial.
#' @param block block size.
#' @param rule rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
#' 'FLGI PM' stands for making decision based on posterior mean;
#' 'FLGI PD' stands for making decision based on posterior distribution;
#' 'CFLGI' stands for controlled forward-looking Gittins Index.
#' @param prior_n a vector representing the number of observations assumed in prior distribution, eg: c(1,1) for a two-armed trial.
#' @param prior_mean1 a vector representing mean of observations assumed in prior distributions, eg: c(0,0,0) for a three-armed trial,
#' rep(0,K) can be used to simplify the process. If a negative effect is expected, adjust the mean to a negative value.
#' @param prior_sd1 a vector representing the standard deviation of observations assumed in prior distribution, eg: rep(1,3) for a three-armed trial.
#' @param stopbound the cut-off value for T test statistics, which is simulated under the null hypothesis.
#' @param mean a vector of mean hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the mean
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the mean for the control and
#' a treatment group, respectively.
#' @param sd a vector of standard deviation hypotheses, for example, as c(0.64,0.64) where 0.64 stands for the standard deviation
#' for both groups. Another example is c(0.64,0.4) where 0.64 and 0.4 stand for the standard deviation for the control and
#' a treatment group, respectively.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_flgi_unknown_var} returns an object of class "flgi". An object of class "flgi" is a list containing 
#' final decision based on the T test statistics with 1 stands for selected and 0 stands for not selected, final decision based on 
#' the maximal Gittins Index value at the final stage, T test statistics, the simulated data set and participants accrued for each arm 
#' at the time of termination of that group in one trial. The simulated data set includes 5 columns: participant ID number, enrollment time, 
#' observed time of results, allocated arm, and participants' result.
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @examples
#' #The forward-looking Gittins Index rule with delayed responses follow a normal 
#' #with a mean of 60 days and a standard deviation of 3 days
#' \donttest{
#' sim_flgi_unknown_var(Gittinstype='UNKV',df=0.5,Pats=10,nMax=50000,
#' TimeToOutcome=expression(rnorm( length( vStartTime ),60, 3)),enrollrate=0.9,
#' K=3,noRuns2=100,Tsize=852,block=20,rule='FLGI PM',prior_n=rep(2,3),
#' prior_mean1=rep(9/100,3),stopbound=(-1.9298),prior_sd1=rep(0.006324555,3),
#' mean=c(9.1/100,8.83/100,8.83/100),sd=c(0.009,0.009,0.009),side='lower')}
#' @references
#' \insertRef{Williamson2019}{RARtrials}


sim_flgi_unknown_var<-function(Gittinstype,df,gittins=NULL,Pats,nMax,TimeToOutcome,enrollrate,K,noRuns2,Tsize,block,rule,
                               prior_n,prior_mean1,prior_sd1,stopbound,mean,sd,side){

  if (is.null(gittins)){
    GI_Normal_unknown <- Gittins(Gittinstype,df)
  }else{
    GI_Normal_unknown <-gittins
  }

  index<-matrix(0,nrow=K,1)
  meanhat<-matrix(0,nrow=1,K)
  sigmahat<-matrix(0,nrow=1,K)
  GI_Std<-rep(0,K)
  zs1<-matrix(0,nrow=1,K-1)
  ap<-matrix(0,nrow=1,K-1)

  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:Tsize], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)

  data1<-matrix(NA_real_,nrow=Tsize,ncol=5)
  data1[,1]<-1:Tsize
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  prior_mean<-prior_mean1
  prior_sd<-prior_sd1
  prior_nn<-rep(0,K)

  prior_n1=prior_n
  nn<-rep(0,K)
  GI<-rep(NA,K)

  for (t in 0:((Tsize/block)-1)){

    alp=allocation_probabilities_unk_var(GI_Normal_unknown=GI_Normal_unknown,tt=t,data1=data1,arms=K,b=block,runs=noRuns2,
                                 posteriormean=prior_mean,posteriornn=prior_nn,posteriorsd=prior_sd,prior_n=prior_n,
                                 prior_mean1=prior_mean1,prior_sd1=prior_sd1,side=side)

    if (rule=='Controlled FLGI'  ){
      alp[1]=1/(K-1)
      elp_e=allocation_probabilities_unk_var1(GI_Normal_unknown=GI_Normal_unknown,tt=t,data1=data1,arms=K,b=block,runs=noRuns2,
                                      posteriormean=prior_mean,posteriornn=prior_nn,posteriorsd=prior_sd,prior_n=prior_n,
                                      prior_mean1=prior_mean1,prior_sd1=prior_sd1,side=side)
      c=alp[1]+sum(elp_e)
      alp=(1/c)*c(alp[1],elp_e)
    }

    alp=cumsum(c(0,alp))
    Pob<-rep(0,block)
    Pos<-rep(0,block)


    for (p in 1:block){

      prior_mean<-prior_mean1
      prior_sd<-prior_sd1

      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          Pos[p]=rnorm(1,mean[k],sd[k])
          data1[t*block+p,4]=k
          data1[t*block+p,5]=Pos[p]
        }
      }

      for (k in 1:K){
        dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[t*block+p,2])),],ncol=5)
        nn[k]=nrow(dataa[which(dataa[,4]==k),,drop=F])
        if (nn[k]>0){
          dataaa<-matrix(dataa[order(dataa[,3]),],ncol=5)
          dataa1<-dataaa[dataaa[,4]==k,5]
          posterior_mean<-rep(NA,nn[k])
          posterior_sd<-rep(NA,nn[k])
          posterior_nn<-rep(NA,nn[k])


          for (kk in 0:(nn[k]-1)){
            posterior_mean[kk+1] <- ((prior_n1[k]+kk)*prior_mean[k]+dataa1[kk+1])/(prior_n1[k]+kk+1)
            posterior_sd[kk+1] <- sqrt(((prior_sd[k])^2)*(prior_n1[k]+kk-1)/(prior_n1[k]+kk) +
                                         (dataa1[kk+1]-prior_mean[k])^2/(prior_n1[k]+kk+1))
            posterior_nn[kk+1]<-prior_n1[k]+kk
            prior_mean[k]<-posterior_mean[kk+1]
            prior_sd[k]<-posterior_sd[kk+1]
          }
          prior_mean[k]<-posterior_mean[nn[k]]
          prior_nn[k]<-posterior_nn[nn[k]]-prior_n1[k]+1
          prior_sd[k]<-posterior_sd[nn[k]]

        }else if (nn[k]==0){
          prior_mean[k]<-prior_mean1[k]
          prior_nn[k]<-0
          prior_sd[k]<-prior_sd1[k]

        }
      }

    }
  }


  if ((Tsize %% block)!=0){
    Pob<-rep(0,(Tsize %% block))
    Posi<-rep(0,(Tsize %% block))


    for (p in 1:(Tsize %% block)){
      prior_mean<-prior_mean1
      prior_sd<-prior_sd1

      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          Posi[p]=rnorm(1,mean=mean[k],sd=sd[k])
          data1[floor(Tsize/block)*block+p,4]=k
          data1[floor(Tsize/block)*block+p,5]=Posi[p]
        }
      }

      for (k in 1:K){
        dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[floor(Tsize/block)*block+p,2])),],ncol=5)
        nn[k]=nrow(dataa[which(dataa[,4]==k),,drop=F])


        if (nn[k]>0){
          dataaa<-matrix(dataa[order(dataa[,3]),],ncol=5)
          dataa1<-dataaa[dataaa[,4]==k ,5]

          posterior_mean<-rep(NA,nn[k])
          posterior_sd<-rep(NA,nn[k])
          posterior_nn<-rep(NA,nn[k])
          for (kk in 0:(nn[k]-1)){
            posterior_mean[kk+1] <- ((prior_n1[k]+kk)*prior_mean[k]+dataa1[kk+1])/(prior_n1[k]+kk+1)
            posterior_sd[kk+1] <- sqrt(((prior_sd[k])^2)*(prior_n1[k]+kk-1)/(prior_n1[k]+kk) +
                                         (dataa1[kk+1]-prior_mean[k])^2/(prior_n1[k]+kk+1))
            posterior_nn[kk+1]<-prior_n1[k]+kk+1
            prior_mean[k]<-posterior_mean[kk+1]
            prior_sd[k]<-posterior_sd[kk+1]
          }
          prior_mean[k]<-posterior_mean[nn[k]]
          prior_nn[k]<-posterior_nn[nn[k]]-prior_n1[k]+1
          prior_sd[k]<-posterior_sd[nn[k]]
        }else if (nn[k]==0){
          prior_mean[k]<-prior_mean1[k]
          prior_nn[k]<-0
          prior_sd[k]<-prior_sd1[k]
        }

      }

    }
  }

  meanhat<-rep(0,K)
  sdhat<-rep(0,K)

  for (k in 1:K){
    nn[k]=nrow(data1[which(data1[,4]==k),,drop=F])
    dataaa<-matrix(data1[order(data1[,3]),],ncol=5)
    dataa1<-dataaa[dataaa[,4]==k,5]

    GI_Std  <- GI_Normal_unknown[ prior_nn[k]+prior_n1[k] ]
    if (side=='upper'){
      GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
    }else if (side=='lower'){
      GI[k]   <- (-prior_mean[k]) + prior_sd[k]*GI_Std
    }

    meanhat[k]=mean(dataa1)
    sdhat[k]=sd(dataa1)
  }


  pc<-matrix(0,1,K-1)

  for (k in 1:(K-1)){
    zs1[1,k]=(meanhat[k+1]-meanhat[1])/sqrt((sdhat[1])^2/nn[1]+(sdhat[k+1])^2/nn[k+1])
    if (is.na(zs1[1,k])){
      zs1[1,k]=0
    }
  }



  for (k in 1:(K-1)){
    if (side=='upper'){
      if(zs1[1,k]>=stopbound ){
        ap[1,k]=1
      }else{
        ap[1,k]=0
      }
    }else if (side=='lower'){
      if(zs1[1,k]<=stopbound){
        ap[1,k]=1
      }else{
        ap[1,k]=0
      }
    }

  }


  #### If options=1, then selection is done based on GI criteria
  indexa<-matrix(0,1,K)

  for (k in 1:K){
    if (rule %in% c('FLGI PM', 'Controlled Gittins','FLGI PD')){
      if (side=='upper'){
        indexa[1,k] =prior_mean[k] + prior_sd[k]*GI[k]
      }else if (side=='lower'){
        indexa[1,k] =(-prior_mean[k]) + prior_sd[k]*GI[k]
      }
    }
  }

  decision= max.col(indexa)

  #return(list(ap,decision,zs1,data1,nn))
  output1<-list(ap,decision,zs1,data1,nn)
  class(output1)<-'flgi'
  
  
  return(output1)
}


#' @export
print.flgi<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[3]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[5]],sep=', ',collapse=', '))
  invisible(x)
}
