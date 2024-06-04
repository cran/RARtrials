#' @title Select au in Bayesian Response-Adaptive Randomization with a Control Group for Binary Endpoint
#' @description \code{brar_select_au_binary} involves selecting au in Bayesian Response-Adaptive Randomization with a control group
#' for binary outcomes with two to five arms. The conjugate prior distributions follow Beta (\eqn{Beta(\alpha,\beta)}) distributions
#' and can be specified individually for each arm.
#' @details This function generates a data set or a value in one iteration for selecting the appropriate au using Bayesian
#' response-adaptive randomization with a control group under null hypotheses with no delay and delayed scenarios.
#' The function can handle trials with up to 5 arms for binary outcomes. This function uses the formula
#' \eqn{\frac{Pr(p_k=max\{p_1,...,p_K\})^{tp}} {\sum_{k=1}^{K}{Pr(p_k=max\{p_1,...,p_K\})^{tp}}}} with \code{side} equals to 'upper',
#' and \eqn{\frac{Pr(p_k=min\{p_1,...,p_K\})^{tp}} {\sum_{k=1}^{K}{Pr(p_k=min\{p_1,...,p_K\}){tp}}}} 
#' with \code{side} equals to 'lower', utilizing available data at each step.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases brar_select_au_binary
#' @export brar_select_au_binary
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
#' @param armn number of total arms in the trial.
#' @param h a vector of success probabilities in hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the success probability
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probabilities
#' for the control and the treatment group, respectively.
#' @param N2 maximal sample size for the trial.
#' @param tp tuning parameter. Some commonly used numbers are 0.5, 1 and n/2N.
#' @param armlabel a vector of arm labels with an example of c(1, 2), where 1 and 2 describe
#' how each arm is labeled in a two-armed trial.
#' @param blocksize size of block used for equal randomization regarding participants in the 'initialization' period.
#' Recommend to be an even multiple of the number of total arms.
#' @param alpha1,beta1 \eqn{\alpha} and \eqn{\beta} in the \eqn{Beta(\alpha,\beta)}, prior for arm 1 which
#' stands for the control. Default value is set to 1.
#' @param alpha2,beta2 \eqn{\alpha} and \eqn{\beta} in the \eqn{Beta(\alpha,\beta)}, prior for arm 2.
#' Default value is set to \code{alpha1} and \code{beta1}.
#' @param alpha3,beta3 \eqn{\alpha} and \eqn{\beta} in the \eqn{Beta(\alpha,\beta)} prior for arm 3.
#' Default value is set to \code{alpha1} and \code{beta1}.
#' @param alpha4,beta4 \eqn{\alpha} and \eqn{\beta} in the \eqn{Beta(\alpha,\beta)} prior for arm 4.
#' Default value is set to \code{alpha1} and \code{beta1}..
#' @param alpha5,beta5 \eqn{\alpha} and \eqn{\beta} in the \eqn{Beta(\alpha,\beta)} prior for arm 5.
#' Default value is set to \code{alpha1} and \code{beta1}.
#' @param minstart a specified number of participants when one starts to check decision rules.
#' @param deltaa a vector of minimal effect expected to be observed for early futility stopping in
#' each arm is approximately \eqn{1\%}. The length of this parameter is \code{armn}-1.
#' @param tpp indicator of \code{tp} equals to n/2N. When \code{tp} is n/2N, \code{tpp} should be assigned 1. Default value is set to 0.
#' @param deltaa1 a vector of pre-specified minimal effect size expected to be observed at the final stage
#' for each arm. The length of this parameter is \code{armn}-1.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param output control the output of \code{brar_select_au_binary}. If user does not specify anything, the function returns
#' the entire data set used to select the stopping boundary for each iteration. If the user specifies 'B', the function
#' only returns the selected stopping boundary for each iteration.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return A list of results generated from formula \eqn{Pr(p_k>p_{control}+\delta|data_{t-1})} at each step.
#' Note that before final stage of the trial, test statistics is calculated from \code{deltaa}, and test statistics is
#' calculated from \code{deltaa1} at the final stage.
#' @importFrom Rdpack reprompt
#' @examples
#' #brar_select_au_binary with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days, where h1=c(0.2,0.4), tp=0.5,
#' #early futility stopping is set at -0.085, and the minimal effect size is 0.1.
#' set.seed(123)
#' stopbound1<-lapply(1:10,function(x){brar_select_au_binary(Pats=10,
#' nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
#' enrollrate=0.1,N1=24,armn=2,h=c(0.3,0.3),N2=224,tp=0.5,armlabel=c(1,2),
#' blocksize=4,alpha1=1,beta1=1,alpha2=1,beta2=1,minstart=24,deltaa=-0.01,
#' tpp=0,deltaa1=0.1,side='upper')})
#' simf<-list()
#' for (xx in 1:10){
#'     if (any(stopbound1[[xx]][24:223,2]<0.01)){
#'       simf[[xx]]<-NA
#'    }  else{
#'       simf[[xx]]<-stopbound1[[xx]][224,2]
#'  }
#'}
#'simf1<-do.call(rbind,simf)
#'sum(is.na(simf1))/10  #1, achieve around 10% futility
#'sum(simf1>0.36,na.rm=TRUE)/10  #0.1
#'#the selected stopping boundary is 0.36 with an overall upper one-sided type I
#'#error of 0.1, based on 10 simulations. It is recommended to conduct more simulations 
#'#(i.e., 1000) to obtain an accurate au.
#' @references 
#' \insertRef{Wathen2017}{RARtrials}


brar_select_au_binary<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,armn,h,N2,tp,armlabel,blocksize,
                          alpha1=1,beta1=1,alpha2=alpha1,beta2=beta1,alpha3=alpha1,beta3=beta1,
                          alpha4=alpha1,beta4=beta1,alpha5=alpha1,beta5=beta1,minstart,deltaa,tpp=0,deltaa1,side,output=NULL,...){

  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)
  assign1<-blockrand(blocksize,N1,armn,armlabel)

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

  armleft<-c(1:armn)
  decision<-rep(NA,armn )
  phi<-rep(NA,armn )
  stopp<-rep(NA,armn )
  simout<-matrix(NA,nrow=N2,ncol=armn)
 
  alpha<-list(alpha1,alpha2,alpha3,alpha4,alpha5)
  beta<-list(beta1,beta2,beta3,beta4,beta5)
  
  for (jjj in minstart:N2){

    if (jjj>minstart){
      treat<-sample(armleft,size =1, prob = as.vector(pii))
      data1[jjj,4]<-treat
      data1[jjj,5]<-rbinom(1,size=1,prob=h[treat])
    }

    if (jjj<N2){
      total<-sum (as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2]))
    }else if (jjj==N2){
      total<-N2
    }

    result<-vector("list",length(armleft))
    mat<-vector("list",armn)

    for (j in 1:length(armleft)) {
      if (total>0){
        if (jjj!=N2){
          data2<-matrix(data1[which(as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2])),],ncol=5)
        }else if (jjj==N2){
          data2<-data1
        }
        tot<-as.numeric(data2[which(data2[,4]==as.numeric(armleft[j])),5])
        mat[[armleft[j]]]<-matrix(c(sum(tot),length(tot) -sum(tot),length(tot)),nrow=1)
      }

    }

    if (length(armleft)>1){
      for (j in 1:length(armleft)){
        if (total>0){
          if (j>1){
            result[[j]]<-pgreater_beta(a1=mat[[1]][1,1]+alpha[[1]],
                                  b1=mat[[1]][1,2]+beta[[1]],
                                  a2=mat[[armleft[j]]][1,1]+alpha[[armleft[j]]],
                                  b2=mat[[armleft[j]]][1,2]+beta[[armleft[j]]],
                                  delta=deltaa[armleft[j]-1],side=side)
          }else if (j==1){
            result[[1]]<-0
          }

        }else if (total==0 ){
          if (j>1){
            result[[j]]<-pgreater_beta(a1=alpha[[1]],
                                  b1=beta[[1]],
                                  a2=alpha[[armleft[j]]],
                                  b2=beta[[armleft[j]]],
                                  delta=deltaa[armleft[j]-1],side=side)
          }else if (j==1){
            result[[1]]<-0
          }
        }

      }

      aloo<-vector("list",length(armleft))
      aloo<-alofun(alpha1=alpha1,beta1=beta1,alpha2=alpha2,beta2=beta2,
                   alpha3=alpha3,beta3=beta3,alpha4=alpha4,beta4=beta4,
                   alpha5=alpha5,beta5=beta5,mat=mat,total=total,armleft=armleft,side=side)

    }


    if (jjj==N2){
      resultt<-vector("list",length(armleft))
      if (length(armleft)>1){
        for (j in 1:length(armleft)){
          if (total>0){
            if (j>1){
              resultt[[j]]<-pgreater_beta(a1=mat[[1]][1,1]+alpha[[1]],
                                     b1=mat[[1]][1,2]+beta[[1]],
                                     a2=mat[[armleft[j]]][1,1]+alpha[[armleft[j]]],
                                     b2=mat[[armleft[j]]][1,2]+beta[[armleft[j]]],
                                     delta=deltaa1[armleft[j]-1],side=side)
            }else if (j==1){
              resultt[[1]]<-0
            }

          }else if (total==0 ){
            if (j>1){
              resultt[[j]]<-pgreater_beta(a1=alpha[[1]],
                                     b1=beta[[1]],
                                     a2=alpha[[armleft[j]]],
                                     b2=beta[[armleft[j]]],
                                     delta=deltaa1[armleft[j]-1],side=side)
            }else if (j==1){
              resultt[[1]]<-0
            }
          }

        }
      }
    }

    pii<-as.data.frame(do.call(cbind,aloo))
    colnames(pii)<-armleft

    if (jjj<N2){
      sim111<-do.call(cbind,result)
      colnames(sim111)<-armleft
      simout[jjj,]<-sim111
    } else if (jjj==N2) {
      sim111t<-do.call(cbind,resultt)
      colnames(sim111t)<-armleft
      simout[jjj,]<-sim111t
    }


    if ( jjj<N2){

      for (yy in 1:length(armleft)){
        if ( pii[colnames(pii) %in%  armleft[yy]]<0.1){
          pii[colnames(pii) %in%  armleft[yy]]=0.1
        }else if ( pii[colnames(pii) %in%  armleft[yy]]>0.9){
          pii[colnames(pii) %in%  armleft[yy]]=0.9
        }
      }


      if (tpp==1){
        pii<-(pii^(jjj/(2*N2)))^tp
      }else if (tpp==0){
        pii<-pii^tp
      }
      pii<-pii/sum(pii)
      pii<-pii[ , order(names(pii))]
      armleft<-sort(armleft,decreasing = FALSE)

    }else if (jjj==N2){

      if (is.null(output)){
        return(simout)
      } else {
        simout1<- rep(NA,armn-1 )
        for (k in 1:(armn-1)){
          if (any(simout[N1:(N2-1),k+1]<0.01)){
            simout1[k]<-NA
          }else{
            simout1[k]<-simout[N2,k+1]
          }
        }
        return(simout1)
      }
    }
  }
}


