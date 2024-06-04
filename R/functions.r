#' @importFrom stats qpois
#' @importFrom stats rexp
#' @importFrom stats runif
#' @importFrom stats rbeta
#' @importFrom stats rbinom
#' @importFrom stats rnorm


###Functions of SimulateAMonthOfAccrualTimes, SimulateArrivalTimes and SimulateOutcomeObservedTime
###are adapted from http://github.com/kwathen/IntroBayesianSimulation with minor modification
SimulateAMonthOfAccrualTimes <- function( dPatsPerMonth , dStartMonth )
{
  nQtyPats    <- 1.2 *qpois(0.9999,dPatsPerMonth)
  vTimes      <- cumsum( rexp( nQtyPats, dPatsPerMonth ) )
  vTimes      <- vTimes[ vTimes < 1 ]
  vTimes      <- vTimes + dStartMonth
  return( vTimes )
}


SimulateArrivalTimes <- function( Pats, nMax )
{
  vTimes <- c()
  if( length( Pats ) == 1 )
  {
    vTimes <- cumsum(rexp(nMax ,Pats))
  }
  else
  {
    dStartMonth <- 0
    nMonth     <- 1
    while( length( vTimes ) < nMax  )
    {
      vTimes      <- c( vTimes, SimulateAMonthOfAccrualTimes( Pats[ nMonth ], dStartMonth ))
      dStartMonth <- dStartMonth + 1

      if( nMonth < length( Pats ) )
        nMonth <- nMonth +  1
    }
    vTimes <- vTimes[ 1:nMax ]
  }
  return( vTimes )
}

SimulateOutcomeObservedTime <- function( vStartTime,TimeToOutcome )
{
  vTimeToOutcome <-eval(TimeToOutcome)
  vObsTime <- vStartTime  + vTimeToOutcome
  return( vObsTime )
}

#### pop simulates the population of patients, returning the start time
#### of becoming ill(initiating study treatment) for all patients, the maximal accrued
#### number of patients in the population, along with the start time of participants enrolled in the trial.
#### Pats: the number of patients accrued within a certain time frame indicates the
#### count of individuals who have been affected by the disease during that specific period,
#### for example, a month or a day. If this number is 10, it represents that
#### 10 people have got the disease within the specified time frame.
#### nMax: the maximal accrued number of patients with the disease, this number
#### should be chosen carefully to ensure a sufficient number of patients are simulated,
#### especially when considering the delay mechanism.
#### enrollrate: probability that patients in the population can enroll in the trial.
#### This parameter is related to the number of people who have been affected by the
#### disease in the population, following an exponential distribution.
pop<-function(Pats,nMax,enrollrate){
  populationtotal<-SimulateArrivalTimes (Pats, nMax)
  vStartTime1<-rbinom(nMax,size=1,enrollrate)
  vStartTime2<- cbind(vStartTime1,populationtotal)
  vStartTime3<-vStartTime2[vStartTime2[,1]==1,]
  return(list(populationtotal,length(populationtotal),as.vector(vStartTime3[,2])))
}

#### blockrand simulates the block equal randomization, mainly used for participants
#### in the burn-in period for response-adaptive randomization.
#### blocksize: size of block used for equal randomization. Recommend to be an even multiple of the
#### number of treatment groups.
#### N: number of patients with equal randomization in the burn-in period.
#### armn: number of total treatment groups in the trial.
#### armlabel: a vector of treatment labels
blockrand = function(blocksize,N,armn,armlabel){
  block = rep(1:ceiling(N/blocksize), each = blocksize)
  a1 = as.data.frame(cbind(block, rand=runif(length(block)), envelope= 1: length(block)))
  a2 = a1[order(a1$block,a1$rand),]
  a2$arm = rep(armlabel,times = length(block)/armn)
  assign = a2[order(a2$envelope),]
  return(assign[,c("block","arm")])
}

#### startfun simulates the equal randomization period before the response-adaptive randomization.
#### It is only used for doubly adaptive biased coin design with minimal variance
#### strategy and maximal power strategy for binary outcomes.
#### popdat: start time of patients enrolled in the trial from the return of
#### third element of pop function.
#### TimeToOutcome Representation of the time distribution of delayed responses. The accrual times
#### could be a month, a week or any other time frame. When the unit changes,
#### the number of TimeToOutcome should also change. It can be in the format
#### of expression(rnorm( length( vStartTime ),30, 3)), representing delayed responses
#### with a normal distribution, where the mean is 30 days and the standard deviation is 3 days.
#### These related functions are adapted from \url{http://github.com/kwathen/IntroBayesianSimulation}.
#### Refer to the website for more details.
#### blocksize: block size.
#### N1: number of participants with equal randomization in the burn-in period.
#### Recommend using 10 percent of the total sample size.
#### armn: number of total treatment groups in the trial.
#### armlabel: armlabel a vector of treatment labels with an example of c(1, 2), where 1 and 2 describes
#### how each treatment group is labeled in a two-armed trial.
#### N2: maximal sample size for the trial.
#### h: a vector of hypothesis, for example, as c(0.1,0.1) where 0.1 stands for the success rate
#### for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success rates
#### for the control and a treatment group, respectively.
startfun<-function(popdat,TimeToOutcome,blocksize,N1,armn,armlabel,N2,h){

  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome )

  assign1<-blockrand(blocksize=blocksize,N=N1,armn=armn,armlabel=armlabel)
  data1<-matrix(NA_real_,nrow=N2,ncol=5)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  data1[1:N1,4]<-assign1$arm[1:N1]

  for (i in 1:(N1)){
    for (j in 1:armn) {
      if (data1[i, 4]==j ){
        data1[i,5]<-rbinom(1,size=1,prob=h[[j]])
      }
    }
  }

  return(data1)
}

#### startfun1 simulates equal randomization period before the response-adaptive randomization.
#### It is only used for A optimal, Aa optimal and analogy RSIHR optimal for continuous outcomes.
#### popdat: output of pop function. It stands for the start time (time of disease) of participants enrolled in the trial.
#### TimeToOutcome Representation of the time distribution of delayed responses. The accrual times
#### could be a month, a week or any other time frame. When the unit changes,
#### the number of TimeToOutcome should also change. It can be in the format
#### of expression(rnorm( length( vStartTime ),30, 3)), representing delayed responses
#### with a normal distribution, where the mean is 30 days and the standard deviation is 3 days.
#### These related functions are adapted from \url{http://github.com/kwathen/IntroBayesianSimulation}.
#### Refer to the website for more details.
#### blocksize: block size.
#### N1: number of participants with equal randomization in the burn-in period.
#### Recommend using 10 percent of the total sample size.
#### armn: number of total treatment groups in the trial.
#### armlabel: armlabel a vector of treatment labels with an example of c(1, 2), where 1 and 2 describes
#### how each treatment group is labeled in a two-armed trial.
#### N2: maximal sample size for the trial.
#### mean: a vector of hypotheses of mean for all treatment groups in the trial, with
#### the first one serving as the control group.
#### sd: a vector of hypotheses of standard deviation for all treatment groups in the trial,
#### with the first one serving as the control group.
startfun1<-function(popdat,TimeToOutcome,blocksize,N1,armn,armlabel,N2,mean,sd){

  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)

  assign1<-blockrand(blocksize=blocksize,N=N1,armn=armn,armlabel=armlabel)
  data1<-matrix(NA_real_,nrow=N2,ncol=5)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  data1[1:N1,4]<-assign1$arm[1:N1]

  for (i in 1:(N1)){
    for (j in 1:armn) {
      if (data1[i, 4]==j ){
        data1[i,5]<-rnorm(1,mean[j],sd[j])
      }
    }
  }

  return(data1)
}

#### alofun calculates the allocation probability of each treatment group in the trial.
#### This function returns a vector of posterior probabilities that treatment group k is maximal
#### considering treatment groups left in the trial. It is only used for Bayesian
#### response-adaptive randomization with a control group using Thall & Wathen method for binary outcomes.
#### alpha1-alpha5: alpha in beta(\alpha,\beta), prior for treatment groups.
#### beta1-beta5: beta in beta(\alpha,\beta), prior for treatment groups.
#### mat: an intermediate variable calculated from brarfun.
#### total: an intermediate variable representing the number of participants with available results.
#### armleft: an intermediate variable representing treatment groups left in the trial at the current stage.
#### side: direction of a one-sided test, with values 'upper' or 'lower'.
alofun<-function(alpha1,beta1,alpha2,beta2,alpha3,beta3,
                 alpha4,beta4,alpha5,beta5,mat,total,armleft,side){
  aloo<-vector("list",length(armleft))
  if (total>0){
    if (length(armleft)==2){
      aloo[[1]]<-pmax_beta(armn=length(armleft),a1=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b1=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a2=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b2=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),side=side)
      aloo[[2]]<-1-aloo[[1]]

    }else if (length(armleft)==3){

      aloo[[1]]<-pmax_beta(armn=length(armleft),a1=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b1=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a2=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b2=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),side=side)

      aloo[[2]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a1=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b1=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),side=side)

      aloo[[3]]<-1-aloo[[1]]-aloo[[2]]

    }else if (length(armleft)==4){

      aloo[[1]]<-pmax_beta(armn=length(armleft),a1=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b1=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a2=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b2=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),side=side)

      aloo[[2]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a1=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b1=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),side=side)

      aloo[[3]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a3=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b3=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a1=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b1=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),side=side)

      aloo[[4]]<-1- aloo[[1]]-aloo[[2]]-aloo[[3]]

    }else if (length(armleft)==5){

      aloo[[1]]<-pmax_beta(armn=length(armleft),a1=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b1=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a2=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b2=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),
                       a5=mat[[armleft[5]]][1,1]+get(sprintf("alpha%s",armleft[5])),
                       b5=mat[[armleft[5]]][1,2]+get(sprintf("beta%s",armleft[5])),side=side)

      aloo[[2]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a1=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b1=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a3=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b3=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),
                       a5=mat[[armleft[5]]][1,1]+get(sprintf("alpha%s",armleft[5])),
                       b5=mat[[armleft[5]]][1,2]+get(sprintf("beta%s",armleft[5])),side=side)

      aloo[[3]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a3=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b3=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a1=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b1=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a4=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b4=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),
                       a5=mat[[armleft[5]]][1,1]+get(sprintf("alpha%s",armleft[5])),
                       b5=mat[[armleft[5]]][1,2]+get(sprintf("beta%s",armleft[5])),side=side)

      aloo[[4]]<-pmax_beta(armn=length(armleft),a2=mat[[1]][1,1]+get(sprintf("alpha%s",1)),
                       b2=mat[[1]][1,2]+get(sprintf("beta%s",1)),
                       a3=mat[[armleft[2]]][1,1]+get(sprintf("alpha%s",armleft[2])),
                       b3=mat[[armleft[2]]][1,2]+get(sprintf("beta%s",armleft[2])),
                       a4=mat[[armleft[3]]][1,1]+get(sprintf("alpha%s",armleft[3])),
                       b4=mat[[armleft[3]]][1,2]+get(sprintf("beta%s",armleft[3])),
                       a1=mat[[armleft[4]]][1,1]+get(sprintf("alpha%s",armleft[4])),
                       b1=mat[[armleft[4]]][1,2]+get(sprintf("beta%s",armleft[4])),
                       a5=mat[[armleft[5]]][1,1]+get(sprintf("alpha%s",armleft[5])),
                       b5=mat[[armleft[5]]][1,2]+get(sprintf("beta%s",armleft[5])),side=side)

      aloo[[5]]<-1-aloo[[1]]-aloo[[2]]-aloo[[3]]-aloo[[4]]
    }
  }else if (total==0 ){
    for(j in 1:length(aloo)){
      aloo[[j]]<-1/length(armleft)
    }
  }
  return(aloo)
}

#### pnichisq_mu calculates the cumulative density function (CDF) of the marginal posterior of mu,
#### which follows a t-distribution from a normal-inverse-gamma distribution. This function is only
#### used in 'pgreater_NIG.r'.
#### x: vector of quantiles.
#### par: vector of current parameters from a normal-inverse-chi-squared distribution.
#### side: direction of a one-sided test, with values 'upper' or 'lower'.
pnichisq_mu<-function(x, par,side){
  m<-(x-par$mu)/sqrt(par$sigsq/par$kappa)
  if (side=='lower'){
    pt(m,df=par$nu,lower.tail=TRUE)  #control>trt
  } else if (side=='upper'){
    pt(m,df=par$nu,lower.tail=FALSE) #control<trt
  }
}

#### pnichisq_mu1 calculates the cumulative density function (CDF) of the marginal posterior of mu,
#### which follows a t-distribution from a normal-inverse-gamma distribution. This function is
#### used in 'pmax.r'.
#### x: vector of quantiles.
#### par: vector of current parameters from a normal-inverse-chi-squared distribution.
#### side: direction of a one-sided test, with values 'upper' or 'lower'.
pnichisq_mu1<-function(x, par,side){
  m<-(x-par$mu)/sqrt(par$sigsq/par$kappa)
  if (side=='lower'){
    pt(m,df=par$nu,lower.tail=FALSE)
  } else if (side=='upper'){
    pt(m,df=par$nu,lower.tail=TRUE)
  }
}

#### dnichisq_mu calculates the cumulative probability density function (PDF) of the marginal posterior of mu,
#### which follows a t-distribution from a normal-inverse-gamma distribution.
#### x: vector of quantiles.
#### par: vector of current parameters from a normal-inverse-chi-squared distribution.
dnichisq_mu<-function(x, par,...){
  m<-(x-par$mu)/sqrt(par$sigsq/par$kappa)
  dt(m,df=par$nu,...)/sqrt(par$sigsq/par$kappa)
}

#### alofun_kn_var calculates the allocation probability of each treatment group in the trial.
#### This function returns a vector of posterior probabilities that treatment group k is maximal
#### considering treatment groups left in the trial. It is only used for Bayesian
#### response-adaptive randomization with a control group using Thall & Wathen method
#### for continuous outcomes with known variance.
#### mat: an intermediate variable calculated from brarfun.
#### total: an intermediate variable representing the number of participants with available results.
#### armleft: an intermediate variable representing treatment groups left in the trial at the current stage.
#### side: direction of a one-sided test, with values 'upper' or 'lower'.
alofun_kn_var<-function(mat,total,armleft,side){
  aloo<-vector("list",length(armleft))
  if (total>0){
    if (length(armleft)==2){
      aloo[[1]]<-pmax_normal(armn=length(armleft),mean1=mat[[1]][1,1],
                                 sd1=mat[[1]][1,3],
                                 mean2=mat[[armleft[2]]][1,1],
                                 sd2=mat[[armleft[2]]][1,3],side=side)
      aloo[[2]]<-1-aloo[[1]]

    }else if (length(armleft)==3){

      aloo[[1]]<-pmax_normal(armn=length(armleft),mean1=mat[[1]][1,1],
                                 sd1=mat[[1]][1,3],
                                 mean2=mat[[armleft[2]]][1,1],
                                 sd2=mat[[armleft[2]]][1,3],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],side=side)

      aloo[[2]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[2]]][1,1],
                                 sd1=mat[[armleft[2]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,3],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],side=side)

      aloo[[3]]<-1-aloo[[1]]-aloo[[2]]

    }else if (length(armleft)==4){

      aloo[[1]]<-pmax_normal(armn=length(armleft),mean1=mat[[1]][1,1],
                                 sd1=mat[[1]][1,2],
                                 mean2=mat[[armleft[2]]][1,1],
                                 sd2=mat[[armleft[2]]][1,3],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],side=side)

      aloo[[2]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[2]]][1,1],
                                 sd1=mat[[armleft[2]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,2],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],side=side)

      aloo[[3]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[3]]][1,1],
                                 sd1=mat[[armleft[3]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,2],
                                 mean3=mat[[armleft[2]]][1,1],
                                 sd3=mat[[armleft[2]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],side=side)

      aloo[[4]]<-1- aloo[[1]]-aloo[[2]]-aloo[[3]]

    }else if (length(armleft)==5){

      aloo[[1]]<-pmax_normal(armn=length(armleft),mean1=mat[[1]][1,1],
                                 sd1=mat[[1]][1,],
                                 mean2=mat[[armleft[2]]][1,1],
                                 sd2=mat[[armleft[2]]][1,3],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],
                                 mean5=mat[[armleft[5]]][1,1],
                                 sd5=mat[[armleft[5]]][1,3],side=side)

      aloo[[2]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[2]]][1,1],
                                 sd1=mat[[armleft[2]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,2],
                                 mean3=mat[[armleft[3]]][1,1],
                                 sd3=mat[[armleft[3]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],
                                 mean5=mat[[armleft[5]]][1,1],
                                 sd5=mat[[armleft[5]]][1,3],side=side)

      aloo[[3]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[3]]][1,1],
                                 sd1=mat[[armleft[3]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,3],
                                 mean3=mat[[armleft[2]]][1,1],
                                 sd3=mat[[armleft[2]]][1,3],
                                 mean4=mat[[armleft[4]]][1,1],
                                 sd4=mat[[armleft[4]]][1,3],
                                 mean5=mat[[armleft[5]]][1,1],
                                 sd5=mat[[armleft[5]]][1,3],side=side)

      aloo[[4]]<-pmax_normal(armn=length(armleft),mean1=mat[[armleft[4]]][1,1],
                                 sd1=mat[[armleft[4]]][1,3],
                                 mean2=mat[[1]][1,1],
                                 sd2=mat[[1]][1,3],
                                 mean3=mat[[armleft[2]]][1,1],
                                 sd3=mat[[armleft[2]]][1,3],
                                 mean4=mat[[armleft[3]]][1,1],
                                 sd4=mat[[armleft[3]]][1,3],
                                 mean5=mat[[armleft[5]]][1,1],
                                 sd5=mat[[armleft[5]]][1,3],side=side)

      aloo[[5]]<-1-aloo[[1]]-aloo[[2]]-aloo[[3]]-aloo[[4]]
    }
  }else if (total==0 ){
    for(j in 1:length(aloo)){
      aloo[[j]]<-1/length(armleft)
    }
  }
  return(aloo)
}

#### alofun_unkn_var calculates the allocation probability of each treatment group in the trial.
#### This function returns a vector of posterior probabilities that treatment group k is maximal
#### considering treatment groups left in the trial. It is only used for Bayesian
#### response-adaptive randomization with a control group using Thall & Wathen method
#### for continuous outcomes with ynknown variance.
#### mat: an intermediate variable calculated from brarfun.
#### total: an intermediate variable representing the number of participants with available results.
#### armleft: an intermediate variable representing treatment groups left in the trial at the current stage.
#### side: direction of a one-sided test, with values 'upper' or 'lower'.
alofun_unk_var<-function(mat,total,armleft,side){
  aloo<-vector("list",length(armleft))
  if (total>0){
    if (length(armleft)==2){
      aloo[[1]]<-pmax_NIX(armn=length(armleft),par1=mat[[1]],par2=mat[[armleft[2]]],side=side)
      aloo[[2]]<-1-aloo[[1]]

    }else if (length(armleft)==3){

      aloo[[1]]<-pmax_NIX(armn=length(armleft),par1=mat[[1]],par2=mat[[armleft[2]]],par3=mat[[armleft[3]]],side=side)

      aloo[[2]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[2]]],par2=mat[[1]],par3=mat[[armleft[3]]],side=side)

      aloo[[3]]<-1-aloo[[1]]-aloo[[2]]

    }else if (length(armleft)==4){

      aloo[[1]]<-pmax_NIX(armn=length(armleft),par1=mat[[1]],par2=mat[[armleft[2]]],par3=mat[[armleft[3]]],par4=mat[[armleft[4]]],side=side)

      aloo[[2]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[2]]],par2=mat[[1]],par3=mat[[armleft[3]]],par4=mat[[armleft[4]]],side=side)

      aloo[[3]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[3]]],par2=mat[[1]],par3=mat[[armleft[2]]],par4=mat[[armleft[4]]],side=side)

      aloo[[4]]<-1- aloo[[1]]-aloo[[2]]-aloo[[3]]

    }else if (length(armleft)==5){

      aloo[[1]]<-pmax_NIX(armn=length(armleft),par1=mat[[1]],par2=mat[[armleft[2]]],par3=mat[[armleft[3]]],par4=mat[[armleft[4]]],par5=mat[[armleft[5]]],side=side)

      aloo[[2]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[2]]],par2=mat[[1]],par3=mat[[armleft[3]]],par4=mat[[armleft[4]]],par5=mat[[armleft[5]]],side=side)

      aloo[[3]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[3]]],par2=mat[[1]],par3=mat[[armleft[2]]],par4=mat[[armleft[4]]],par5=mat[[armleft[5]]],side=side)

      aloo[[4]]<-pmax_NIX(armn=length(armleft),par1=mat[[armleft[4]]],par2=mat[[1]],par3=mat[[armleft[2]]],par4=mat[[armleft[3]]],par5=mat[[armleft[5]]],side=side)

      aloo[[5]]<-1-aloo[[1]]-aloo[[2]]-aloo[[3]]-aloo[[4]]
    }
  }else if (total==0 ){
    for(j in 1:length(aloo)){
      aloo[[j]]<-1/length(armleft)
    }
  }
  return(aloo)
}


#### allocation_probabilities simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block.
#### It is only used for forward-looking Gittins index algorithm with binary outcomes.
#### GI_binary: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### K1: number of total treatment groups in the trial.
#### I0: a matrix with K rows and 2 columns, where the numbers inside are equal to the prior parameters, and
#### K is equal to the total number of treatment groups. For example, matrix(1,nrow=2,ncol=2) means that the prior
#### distributions for two-armed trials are beta(1,1); matrix(c(2,3),nrow=2,ncol=2,byrow = TRUE) means that the prior
#### distributions for two-armed trials are beta(2,3). The first column represents the prior of the number of successes,
#### and the second column represents the prior of the number of failures.
#### block: block size
#### noRuns2: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### rule: rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
#### 'FLGI PM' stands for making decision based on posterior mean;
#### 'FLGI PD' stands for making decision based on posterior distribution;
#### 'CFLGI' stands for controlled forward looking Gittins index.
allocation_probabilities<-function(GI_binary,tt,data1,K1,I0,block,noRuns2,rule){

  index<-matrix(0,nrow=K1,1)
  selected<-matrix(0,nrow=noRuns2,block)
  prob<-matrix(0,K1,block)

  for (j in 1:noRuns2) {
    n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
    s=matrix(I0[,1]+1,nrow=nrow(I0),1)
    f=matrix(I0[,2]+1,nrow=nrow(I0),1)

    for (t1 in 0: (block-1)){
      for (k in 1:K1){
        index[k,1]=GI_binary[f[k,1],s[k,1]]
      }
      max_index=max(index)
      kmax=min(match(max(index), index))
      idx= which(as.vector(index)==max(index))

      if (length(idx)>1){
        posi=sample(1:length(idx),1)
        kmax=idx[posi]
        max_index=as.vector(index)[kmax]
      }
      selected[j,t1+1]=kmax
      snext=s
      fnext=f
      nnext=n

      if (rule=='FLGI PM' |rule=='Controlled FLGI'){
        probSuc_kmax=s[kmax,1]/(s[kmax,1]+f[kmax,1])
        if(runif(1)<=probSuc_kmax){
          Pos=1
        }else{
          Pos=0
        }
      }else if (rule=='FLGI PD'){
        probSuc_kmax=rbeta(1,s[kmax,1],f[kmax,1])
        Pos= rbinom(1,1,probSuc_kmax)
      }
      nnext[kmax,1]=n[kmax,1]+1

      data1[tt*block+t1+1,4]=kmax
      data1[tt*block+t1+1,5]=Pos
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2]))

      for (k in 1:K1){
        if (total1>0){
          dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2])),],ncol=5)
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

  for (i in 1:block){
    for(k in 1:K1){
      prob[k,i]=sum(selected[,i]==k)/noRuns2
    }
  }

  allocation_probabilities=rowMeans(prob)

  return(allocation_probabilities)
}

#### allocation_probabilities1 simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block.
#### It is only used for controlled forward-looking Gittins index algorithm with binary outcomes.
#### GI_binary: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### K1: number of total treatment groups in the trial.
#### I0: a matrix with K rows and 2 columns, where the numbers inside are equal to the prior parameters, and
#### K is equal to the total number of treatment groups. For example, matrix(1,nrow=2,ncol=2) means that the prior
#### distributions for two-armed trials are beta(1,1); matrix(c(2,3),nrow=2,ncol=2,byrow = TRUE) means that the prior
#### distributions for two-armed trials are beta(2,3). The first column represents the prior of the number of successes,
#### and the second column represents the prior of the number of failures.
#### block: block size
#### noRuns2: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### rule: rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
#### 'FLGI PM' stands for making decision based on posterior mean;
#### 'FLGI PD' stands for making decision based on posterior distribution;
#### 'CFLGI' stands for controlled forward looking Gittins index.
allocation_probabilities1<-function(GI_binary,tt,data1,K1,I0,block,noRuns2,rule){

  index<-matrix(0,nrow=K1,1)
  selected<-matrix(0,nrow=noRuns2,block)
  prob<-matrix(0,K1,block)

  for (j in 1:noRuns2) {
    n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
    s=matrix(I0[,1]+1,nrow=nrow(I0),1)
    f=matrix(I0[,2]+1,nrow=nrow(I0),1)

    for (t1 in 0: (block-1)){
      for (k in 1:K1){
        index[k,1]=GI_binary[f[k,1],s[k,1]]
      }
      max_index=max(index)
      kmax=min(match(max(index), index))
      idx= which(as.vector(index)==max(index))
      if (length(idx)>1){
        posi=sample(1:length(idx),1)
        kmax=idx[posi]
        max_index=as.vector(index)[kmax]
      }
      selected[j,t1+1]=kmax
      snext=s
      fnext=f
      nnext=n
      if (rule=='FLGI PM'|rule=='Controlled FLGI'){
        probSuc_kmax=s[kmax,1]/(s[kmax,1]+f[kmax,1])
        if(runif(1)<=probSuc_kmax){
          Pos=1
        }else{
          Pos=0
        }
      }else if (rule=='FLGI PD'){
        probSuc_kmax=rbeta(1,s[kmax,1],f[kmax,1])
        Pos= rbinom(1,1,probSuc_kmax)
      }
      nnext[kmax,1]=n[kmax,1]+1
      data1[tt*block+t1+1,4]=kmax+1
      data1[tt*block+t1+1,5]=Pos
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2]))

      for (k in 1:K1){
        if (total1>0){
          dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2])),],ncol=5)
          snext[k,1]=nrow(dataa[dataa[,4]==k+1 & dataa[,5]==1,,drop=F])+2
          fnext[k,1]=nrow(dataa[dataa[,4]==k+1 & dataa[,5]==0,,drop=F])+2

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

  for (i in 1:block){
    for(k in 1:K1){
      prob[k,i]=sum(selected[,i]==k)/noRuns2
    }
  }

  allocation_probabilities=rowMeans(prob)

  return(allocation_probabilities)
}


#### allocation_probabilities_unk_var simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block.
#### It is only used for forward-looking Gittins index algorithm with continuous outcomes
#### and unknown variances.
#### GI_Normal_unknown: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### arms: number of total treatment groups in the trial.
#### b: block size.
#### runs: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### posteriormean: a vector of the calculated posterior means for each treatment group at the current stage.
#### posteriornn: a vector of the calculated posterior number of participants for each treatment group at the current stage.
#### posteriorsd: a vector of the calculated posterior standard deviations for each treatment group at the current stage.
#### prior_n: a vector representing the number of observations assumed in prior distributions, eg: c(1,1) for a two-armed trial.
#### prior_mean1: a vector representing the mean of observations assumed in prior distributions, eg: c(0,0,0) for a three-armed trial.
#### If a negative effect is expected, adjust the mean to a negative value; rep(0,K) can be used to simplify the process.
#### prior_sd1: a vector representing the standard deviation of observations assumed in prior distributions, eg: rep(1,3) for a three-armed trial.
#### side: direction of one-sided test with the values of 'upper' or 'lower'.
allocation_probabilities_unk_var <- function(GI_Normal_unknown,tt,data1,arms,b,runs,
                                             posteriormean,posteriornn,posteriorsd,prior_n,
                                             prior_mean1,prior_sd1,side){


  action <- matrix(0, runs, arms)
  posterior_sd1=posteriorsd
  posterior_mean1=posteriormean
  posterior_nn1=posteriornn
  prior_n1=prior_n

  for (i in 1:runs){
    GI<-rep(NA,arms)
    nn<-rep(0,arms)

    for (k in 1:arms){
      GI_Std <- GI_Normal_unknown[ posterior_nn1[k] +prior_n1[k]]
      if (side=='upper'){
        GI[k]  <- posterior_mean1[k] + posterior_sd1[k]*GI_Std
      }else if (side=='lower'){
        GI[k]  <- (-posterior_mean1[k]) + posterior_sd1[k]*GI_Std
      }
    }




    for (pts in 1:b){

      optimal_action <- which.is.max(GI)
      action[ i , optimal_action ] <- action[ i , optimal_action ] + 1
      data1[tt*b+pts,4]=optimal_action
      if (pts==1){
        outcome <- rnorm(1,mean=posterior_mean1[optimal_action],sd=posterior_sd1[optimal_action])
      }else{
        outcome <- rnorm(1,mean=prior_mean[optimal_action],sd=prior_sd[optimal_action])
      }


      prior_mean<-prior_mean1
      prior_nn<-rep(0,arms)
      prior_sd<-prior_sd1
      data1[tt*b+pts,5]=outcome
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b+pts,2])),],ncol=5)

      for (k in 1:arms){
        dataaa<-matrix(dataa[order(dataa[,3]),],ncol=5)
        nn[k]=nrow(dataaa[dataaa[,4]==k ,,drop=F])
        if (nn[k]>0){
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
          prior_nn[k]<-posterior_nn[nn[k]]
          prior_sd[k]<-posterior_sd[nn[k]]
          GI_Std  <- GI_Normal_unknown[ prior_nn[k]+1 ]
          if (side=='upper'){
            GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
          }else if (side=='lower'){
            GI[k]   <- (-prior_mean[k]) + prior_sd[k]*GI_Std
          }

        }else if (nn[k]==0){

          prior_mean[k]<-prior_mean1[k]
          prior_nn[k]<-prior_n1[k]
          prior_sd[k]<-prior_sd1[k]
          GI_Std  <- GI_Normal_unknown[prior_nn[k]]
          if (side=='upper'){
            GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
          }else if (side=='lower'){
            GI[k]   <- (-prior_mean[k]) + prior_sd[k]*GI_Std
          }

        }
      }

    }
  }
  flgi_probs <- apply(action, 2, mean) / b
  return(flgi_probs)

}

#### allocation_probabilities_unk_var1 simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block.
#### It is only used for controlled forward-looking Gittins index algorithm with continuous outcomes
#### and unknown variances.
#### GI_Normal_unknown: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### arms: number of total treatment groups in the trial.
#### b: block size.
#### runs: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### posteriormean: a vector of the calculated posterior means for each treatment group at the current stage.
#### posteriornn: a vector of the calculated posterior number of participants for each treatment group at the current stage.
#### posteriorsd: a vector of the calculated posterior standard deviations for each treatment group at the current stage.
#### prior_n: a vector representing the number of observations assumed in prior distributions, eg: c(1,1) for a two-armed trial.
#### prior_mean1: a vector representing the mean of observations assumed in prior distributions, eg: c(0,0,0) for a three-armed trial.
#### If a negative effect is expected, adjust the mean to a negative value; rep(0,K) can be used to simplify the process.
#### prior_sd1: a vector representing the standard deviation of observations assumed in prior distributions, eg: rep(1,3) for a three-armed trial.
#### side: direction of one-sided test with the values of 'upper' or 'lower'.
allocation_probabilities_unk_var1 <- function(GI_Normal_unknown,tt,data1,arms,b,runs,
                                              posteriormean,posteriornn,posteriorsd,prior_n,
                                              prior_mean1,prior_sd1,side){

  action <- matrix(0, runs, arms)
  posterior_sd1=posteriorsd
  posterior_mean1=posteriormean
  posterior_nn1=posteriornn
  prior_n1=prior_n

  for (i in 1:runs){
    GI<-rep(NA,arms)
    nn<-rep(0,arms)

    for (k in 2:arms){
      GI_Std <- GI_Normal_unknown[ posterior_nn1[k] +prior_n1[k]]
      if (side=='upper'){
        GI[k]  <- posterior_mean1[k] + posterior_sd1[k]*GI_Std
      }else if (side=='lower'){
        GI[k]  <- (-posterior_mean1[k]) + posterior_sd1[k]*GI_Std
      }
    }




    for (pts in 1:b){


      optimal_action <- which.is.max(!is.na(GI))+1
      action[ i , optimal_action ] <- action[ i , optimal_action ] + 1
      data1[tt*b+pts,4]=optimal_action
      if (pts==1){
        outcome <- rnorm(1,mean=posterior_mean1[optimal_action],sd=posterior_sd1[optimal_action])
      }else{
        outcome <- rnorm(1,mean=prior_mean[optimal_action],sd=prior_sd[optimal_action])
      }
      prior_mean<-prior_mean1
      prior_nn<-rep(0,arms)
      prior_sd<-prior_sd1
      data1[tt*b+pts,5]=outcome
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b+pts,2])),],ncol=5)

      for (k in 2:arms){
        dataaa<-matrix(dataa[order(dataa[,3]),],ncol=5)
        nn[k]=nrow(dataaa[dataaa[,4]==k ,,drop=F])
        if (nn[k]>0){
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
          prior_nn[k]<-posterior_nn[nn[k]]
          prior_sd[k]<-posterior_sd[nn[k]]
          GI_Std  <- GI_Normal_unknown[ prior_nn[k]+1 ]

          if (side=='upper'){
            GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
          }else if (side=='lower'){
            GI[k]   <- (-prior_mean[k]) + prior_sd[k]*GI_Std
          }

        }else if (nn[k]==0){

          prior_mean[k]<-prior_mean1[k]
          prior_nn[k]<-prior_n1[k]
          prior_sd[k]<-prior_sd1[k]
          GI_Std  <- GI_Normal_unknown[prior_nn[k]]

          if (side=='upper'){
            GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
          }else if (side=='lower'){
            GI[k]   <- (-prior_mean[k]) + prior_sd[k]*GI_Std
          }

        }
      }

    }
  }
  flgi_probs <- apply(action, 2, mean) / b
  return(flgi_probs[2:arms])

}


#### allocation_probabilities_kn_var simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block. It is only used for forward-looking Gittins index
#### algorithm with continuous outcomes and known variances.
#### GI_Normal_known: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### arms: number of total treatment groups in the trial.
#### b: block size.
#### runs: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### prior_mean: a vector representing the mean of observations assumed in prior distribution,
#### eg: c(0,0,0) for a three-armed trial. If a negative effect is expected, adjust the
#### mean to a negative value; rep(0,K) can be used to simplify the process.
#### prior_n: a vector representing the number of observations assumed in prior distribution, eg: c(1,1) for a two-armed trial.
#### sd1: a vector of standard deviation hypothesis, for example, as c(0.64,0.64) where 0.64 stands for the standard deviation
#### for both groups. Another example is c(0.64,0.4) where 0.64 and 0.4 stand for the standard deviation for control and
#### a treatment group, respectively.
#### side: direction of one-sided test with the values of 'upper' or 'lower'.
allocation_probabilities_kn_var <- function(GI_Normal_known,tt,data1,arms,b,runs,prior_mean,prior_n,sd1,side){

  action <- matrix(0, runs, arms)
  nn1<-rep(NA,arms)
  posterior_mean<-rep(NA,arms)
  GI<-rep(NA,arms)
  sample_mean<-rep(NA,arms)

  for (k in 1:arms){
    nn1[k]=nrow(data1[data1[,4]==k,,drop=F])
    if (nn1[k]>0){
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b,2])),],ncol=5)
      nn1[k]=nrow(dataa[dataa[,4]==k,,drop=F])+prior_n[k]
      posterior_mean[k]<-(sum(dataa[dataa[,4]==k,5])+prior_mean[k]*prior_n[k]) / nn1[k]
    }else if (nn1[k]==0){
      nn1[k]=prior_n[k]
      posterior_mean[k]<-prior_mean[k]
    }
  }


  for (i in 1:runs){

    GI <- rep(NA, arms)
    nn<-rep(NA,arms)
    for (k in 1:arms){
     # indexs<-nn[k]
      GI_Std <- GI_Normal_known[ nn1[k]]
      if (side=='upper'){
        GI[k]  <- posterior_mean[k] + sd1[k]*GI_Std
      } else if (side=='lower'){
        GI[k]  <- (-posterior_mean[k]) + sd1[k]*GI_Std
      }
    }


    for (pts in 1:b){

      optimal_action <- which.is.max(GI)
      action[ i , optimal_action ] <- action[ i , optimal_action ] + 1
      data1[tt*b+pts,4]=optimal_action
      if (pts==1){
        outcome <- rnorm(1, posterior_mean[optimal_action], sd1[optimal_action])
      }else{
        outcome <- rnorm(1,sample_mean[optimal_action], sd1[optimal_action])
      }
      data1[tt*b+pts,5]=outcome
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b+pts,2])),],ncol=5)

      for (k in 1:arms){
        nn[k]=nrow(dataa[dataa[,4]==k,,drop=F])
        if (nn[k]>0){
          nn[k]=nrow(dataa[dataa[,4]==k,,drop=F])+prior_n[k]
          sample_mean[k] <- (sum(dataa[dataa[,4]==k,5])+prior_mean[k]*prior_n[k]) / nn[k]
          indexs<-nn[k]
          GI_Std  <- GI_Normal_known[ indexs ]
          if (side=='upper'){
             GI[k]   <- sample_mean[k] + sd1[k]*GI_Std
          } else if (side=='lower'){
            GI[k]   <- (-sample_mean[k]) + sd1[k]*GI_Std
          }
        }else if (nn[k]==0){

          nn[k]=prior_n[k]
          sample_mean[k] <- prior_mean[k]
          GI_Std  <- GI_Normal_known[ nn[k] ]
          if (side=='upper'){
            GI[k]   <- sample_mean[k] + sd1[k]*GI_Std
          } else if (side=='lower'){
            GI[k]   <- (-sample_mean[k]) + sd1[k]*GI_Std
          }
        }
      }
    }
  }
  flgi_probs <- apply(action, 2, mean) / b
  return(flgi_probs)

}

#### allocation_probabilities_kn_var1 simulates allocation probabilities for each participant to different
#### treatment groups in the current block, returning the allocation probabilities for
#### each participant in the current block.
#### It is only used for controlled forward-looking Gittins index algorithm with continuous outcomes
#### and known variances.
#### GI_Normal_known: Gittins index table used for analysis.
#### tt: accumulated number of blocks up to the current block.
#### data1: the simulated data set involving participants in the trial with available
#### information at this block.
#### arms: number of total treatment groups in the trial.
#### b: block size.
#### runs: number of simulations for simulated allocation probabilities within each block. Default value is
#### set to 100 times, which is recommended in Villar et al., 2015.
#### prior_mean: a vector representing the mean of observations assumed in prior distribution,
#### eg: c(0,0,0) for a three-armed trial. If a negative effect is expected, adjust the
#### mean to a negative value; rep(0,K) can be used to simplify the process.
#### prior_n: a vector representing the number of observations assumed in prior distribution, eg: c(1,1) for a two-armed trial.
#### sd1: a vector of standard deviation hypothesis, for example, as c(0.64,0.64) where 0.64 stands for the standard deviation
#### for both groups. Another example is c(0.64,0.4) where 0.64 and 0.4 stand for the standard deviation for control and
#### a treatment group, respectively.
#### side: direction of one-sided test with the values of 'upper' or 'lower'.
allocation_probabilities_kn_var1 <- function(GI_Normal_known,tt,data1,arms,b,runs,prior_mean,prior_n,sd1,side){

  action <- matrix(0, runs, arms)
  nn1<-rep(NA,arms)
  posterior_mean<-rep(NA,arms)
  GI<-rep(NA,arms)
  sample_mean<-rep(NA,arms)

  for (k in 1:arms){
    nn1[k]=nrow(data1[data1[,4]==k,,drop=F])
    if (nn1[k]>0){
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b,2])),],ncol=5)
      nn1[k]=nrow(dataa[dataa[,4]==k,,drop=F])+prior_n[k]
      posterior_mean[k]<-(sum(dataa[dataa[,4]==k,5])+prior_mean[k]*prior_n[k]) / nn1[k]
    }else if (nn1[k]==0){
      nn1[k]=prior_n[k]
      posterior_mean[k]<-prior_mean[k]
    }
  }


  for (i in 1:runs){

    GI <- rep(NA, arms)
    nn<-rep(NA,arms)
    for (k in 2:arms){
     # indexs<-nn[k]
      GI_Std <- GI_Normal_known[ nn1[k]]
      if (side=='upper'){
        GI[k]  <- posterior_mean[k] + sd1[k]*GI_Std
      } else if (side=='lower'){
        GI[k]  <- (-posterior_mean[k]) + sd1[k]*GI_Std
      }
    }


    for (pts in 1:b){

      optimal_action <- which.is.max(GI[!is.na(GI)])+1
      action[ i , optimal_action ] <- action[ i , optimal_action ] + 1
      data1[tt*b+pts,4]=optimal_action
      if (pts==1){
        outcome <- rnorm(1, posterior_mean[optimal_action], sd1[optimal_action])
      }else{
        outcome <- rnorm(1,sample_mean[optimal_action], sd1[optimal_action])
      }
      data1[tt*b+pts,5]=outcome
      dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*b+pts,2])),],ncol=5)

      for (k in 2:arms){
        nn[k]=nrow(dataa[dataa[,4]==k,,drop=F])
        if (nn[k]>0){
          nn[k]=nrow(dataa[dataa[,4]==k,,drop=F])+prior_n[k]
          sample_mean[k] <- (sum(dataa[dataa[,4]==k,5])+prior_mean[k]*prior_n[k]) / nn[k]
          indexs<-nn[k]
          GI_Std  <- GI_Normal_known[ indexs ]
          if (side=='upper'){
            GI[k]   <- sample_mean[k] + sd1[k]*GI_Std
          } else if (side=='lower'){
            GI[k]   <- (-sample_mean[k]) + sd1[k]*GI_Std
          }

        }else if (nn[k]==0){

          nn[k]=prior_n[k]
          sample_mean[k] <- prior_mean[k]
          GI_Std  <- GI_Normal_known[ nn[k] ]
          if (side=='upper'){
            GI[k]   <- sample_mean[k] + sd1[k]*GI_Std
          } else if (side=='lower'){
            GI[k]   <- (-sample_mean[k]) + sd1[k]*GI_Std
          }

        }

      }
    }
  }
  flgi_probs <- apply(action, 2, mean) / b
  return(flgi_probs[2:arms])

}

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


# #### different print methods
# print.aoptimal<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste((round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# 
# print.dabcd<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste((round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# 
# print.brar<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# 
# print.flgi<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste(round(x[[3]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[5]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# 
# print.RSIHRoptimal<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# print.rptw<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
# print.Aaoptimal<-function(x,...){
#   cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
#   cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
#   cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
#   invisible(x)
# }
# 
