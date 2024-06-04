#' @title Simulate a Trial Using Forward-Looking Gittins Index for Binary Endpoint
#' @description Function for simulating a trial using the forward-looking Gittins Index rule and the controlled forward-looking
#' Gittins Index rule for binary outcomes in trials with 2-5 arms. The conjugate prior distributions
#' follow Beta (\eqn{Beta(\alpha,\beta)}) distributions and should be the same for each arm.
#' @details This function simulates a trial using the forward-looking Gittins Index rule or the
#' controlled forward-looking Gittins Index rule under both no delay and delayed scenarios.
#' The cut-off value used for \code{stopbound} is obtained by simulations using \code{flgi_stop_bound_binary}.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_flgi_binary
#' @export sim_flgi_binary
#' @param Gittinstype type of Gittins indices, should be set to 'binary' in this function.
#' @param df discount factor which is the multiplier for loss at each additional patient in the future.
#' Available values are 0, 0.5, 0.7, 0.99 and 0.995. The maximal sample size can be up to 2000.
#' @param gittins user specified Gittins indices for calculation in this function. Recommend using the
#' \code{bmab_gi_multiple_ab} function from \code{gittins} package. If \code{gittins} is provided,
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
#' @param I0 a matrix with K rows and 2 columns, where the numbers inside are equal to the prior parameters, and
#' K is equal to the total number of arms. For example, matrix(1,nrow=2,ncol=2) means that the prior
#' distributions for two-armed trials are beta(1,1); matrix(c(2,3),nrow=2,ncol=2,byrow = TRUE) means that the prior
#' distributions for two-armed trials are beta(2,3). The first column represents the prior of the number of successes,
#' and the second column represents the prior of the number of failures.
#' @param K number of total arms in the trial.
#' @param noRuns2 number of simulations for simulated allocation probabilities within each block. Default value is
#' set to 100, which is recommended in \insertCite{Villar2015}{RARtrials}.
#' @param Tsize maximal sample size for the trial.
#' @param ptrue a vector of hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the success probability
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probability for the control and
#' the treatment group, respectively.
#' @param block block size.
#' @param rule rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
#' 'FLGI PM' stands for making decision based on posterior mean;
#' 'FLGI PD' stands for making decision based on posterior distribution;
#' 'CFLGI' stands for controlled forward-looking Gittins Index.
#' @param ztype Z test statistics, with choice of values from 'pooled' and 'unpooled'.
#' @param stopbound the cut-off value for Z test statistics, which is simulated under the null hypothesis.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_flgi_binary} returns an object of class "flgi". An object of class "flgi" is a list containing 
#' final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected, final decision based on 
#' the maximal Gittins Index value at the final stage, Z test statistics, the simulated data set and participants accrued for each arm 
#' at the time of termination of that group in one trial. The simulated data set includes 5 columns: participant ID number, enrollment time, 
#' observed time of results, allocated arm, and participants' result.
#' @importFrom stats runif
#' @examples
#' #The forward-looking Gittins Index rule with delayed responses follow a normal distribution
#' #with a mean of 60 days and a standard deviation of 3 days
#' \donttest{
#' sim_flgi_binary(Gittinstype='Binary',df=0.5,Pats=10,nMax=50000,TimeToOutcome=expression(
#' rnorm( length( vStartTime ),60, 3)),enrollrate=0.9,I0= matrix(1,nrow=2,2),
#' K=2,Tsize=992,ptrue=c(0.6,0.7),block=20,rule='FLGI PM',ztype='unpooled',
#' stopbound=1.9991,side='upper')}
#' @references 
#' \insertRef{Villar2015}{RARtrials}

sim_flgi_binary<-function(Gittinstype,df,gittins=NULL,Pats,nMax,TimeToOutcome,enrollrate,I0,K,noRuns2=100,Tsize,ptrue,block,rule,ztype,stopbound,side){

  if (is.null(gittins)){
    GI_binary <- Gittins(Gittinstype,df)
  }else{
    GI_binary <- gittins
  }

  index<-matrix(0,nrow=K,1)
  phat<-matrix(0,nrow=1,K)
  sigmahat<-matrix(0,nrow=1,K)
  ns<-matrix(0,nrow=1,K)
  sn<-matrix(0,nrow=1,K)
  zs1<-matrix(0,nrow=1,K-1)
  ap<-matrix(0,nrow=1,K-1)


  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:Tsize], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)

  data1<-matrix(NA_real_,nrow=Tsize,ncol=5)
  data1[,1]<-1:Tsize
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime

  n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
  s=matrix(I0[,1]+1,nrow=nrow(I0),1)
  f=matrix(I0[,2]+1,nrow=nrow(I0),1)

  for (t in 0:((Tsize/block)-1)){

    alp=allocation_probabilities(GI_binary=GI_binary,tt=t,data1=data1,I0=cbind(s-1,f-1),block=block,noRuns2=noRuns2,K1=K,rule=rule)
    if (rule=='Controlled FLGI'  ){
      alp[1]=1/(K-1)
      elp_e=allocation_probabilities1(GI_binary=GI_binary,tt=t,data1=data1,I0=cbind(s[2:K,]-1,f[2:K,]-1),block=block,noRuns2=noRuns2,K1=K-1,rule='FLGI PM')
      c=alp[1]+sum(elp_e)
      alp=(1/c)*c(alp[1],elp_e)
    }

    alp=cumsum(c(0,alp))

    snext=s
    fnext=f
    nnext=n

    Pob<-rep(0,block)
    Pos<-rep(0,block)
    for (p in 1:block){
      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          nnext[k]=n[k]+1
          if (runif(1)<=ptrue[k]){
            Pos[p]=1
          }else{
            Pos[p]=0
          }
          data1[t*block+p,4]=k
          data1[t*block+p,5]=Pos[p]
        }
      }
          total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[t*block+p,2]))

      for (k in 1:K){
          if (total1>0){
            dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[t*block+p,2])),],ncol=5)
            snext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==1,,drop=F])+2
            fnext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==0,,drop=F])+2
          }else if (total1==0){
            snext[k,1]=s[k,1]
            fnext[k,1]=f[k,1]
          }
        }

      s=snext
      f=fnext
      n=nnext
    }
  }


  if ((Tsize %% block)!=0){
    Pob<-rep(0,(Tsize %% block))
    Posi<-rep(0,(Tsize %% block))
    for (p in 1:((Tsize %% block))){
      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          nnext[k]=n[k]+1
          if (runif(1)<=ptrue[k]){
            Posi[p]=1
          }else {
            Posi[p]=0
          }
          data1[floor(Tsize/block)*block+p,4]=k
          data1[floor(Tsize/block)*block+p,5]=Posi[p]
        }
      }
          total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[floor(Tsize/block)*block+p,2]))

        for (k in 1:K){
          if (total1>0){
            dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[floor(Tsize/block)*block+p,2])),],ncol=5)
            snext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==1,,drop=F])+2
            fnext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==0,,drop=F])+2
          }else if (total1==0){
            snext[k,1]=s[k,1]
            fnext[k,1]=f[k,1]
          }
        }

      s=snext
      f=fnext
      n=nnext
    }
  }

  for (k in 1:K){
    s[k,1]=nrow(data1[data1[,4]==k & data1[,5]==1,,drop=F])+2
    f[k,1]=nrow(data1[data1[,4]==k & data1[,5]==0,,drop=F])+2
    n[k,1]=nrow(data1[data1[,4]==k ,,drop=F])+4
  }


  ns[1,]=n-2
  sn[1,]=s-1
  phat[1,]=(s-1)/(n-2)

  if (ztype=='unpooled'){
    sigmahat[1,]=(phat[1,]*(1-phat[1,]))/ns[1,]
  } else if (ztype=='pooled'){
    for (k in 2:K){
      sigmahat[1,k]= (sum(sn[1]+sn[k])/sum(ns[1]+ns[k]))*
        (1-(sum(sn[1]+sn[k])/sum(ns[1]+ns[k])))*
        (1/ns[1] +1/ns[k])
    }
  }


  sigma<-matrix(0,K-1,K-1)
  sigmat<-matrix(0,K-1,K-1)
  pc<-matrix(0,1,K-1)

  for (k in 1:(K-1)){
    if (ztype=='unpooled'){
      zs1[1,k]=(phat[1,k+1]-phat[1,1])/sqrt(sigmahat[1,1]+sigmahat[1,k+1])
    } else if (ztype=='pooled'){
      zs1[1,k]=(phat[1,k+1]-phat[1,1])/sqrt(sigmahat[1,k+1])
    }
  }

  b1<-matrix(0,nrow=1,(K-1))

  for (k in 1:(K-1)){
    if (side=='upper'){
      if(zs1[1,k]>=stopbound ){
         b1[1,k]=1
      }else{
         b1[1,k]=0
      }

    }else if (side=='lower'){
      if(zs1[1,k]<=stopbound ){
        b1[1,k]=1
      }else{
        b1[1,k]=0
      }
    }
  }


  indexa<-matrix(0,1,K)

  for (k in 1:K){
    indexa[1,k] = GI_binary[ns[1,k]-sn[1,k]+2,sn[1,k]+1]
  }
  decision=max.col(indexa)

  #return(list(b1,decision,zs1,data1,n[,1]-4))
  output1<-list(b1,decision,zs1,data1,n[,1]-4)
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