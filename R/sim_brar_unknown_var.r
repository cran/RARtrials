#' @title Simulate a Trial Using Bayesian Response-Adaptive Randomization with a Control Group for Continuous Endpoint with Unknown Variances
#' @description \code{sim_brar_unknown_var} simulate a trial with two to five arms using Bayesian Response-Adaptive 
#' Randomization with a control group for continuous outcomes with unknown variances. The conjugate prior distributions
#' follow Normal-Inverse-Gamma (NIG) (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)}) 
#' distributions and can be specified individually for each arm.
#' @details This function generates a designed trial using Bayesian response-adaptive randomization with
#' a control group under no delay and delayed scenarios for continuous outcomes with unknown variances. 
#' The function can handle trials with up to 5 arms. This function uses the formula
#' \eqn{\frac{Pr(\mu_k=max\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k=max\{\mu_1,...,\mu_K\})^{tp}}}} with \code{side} equals to 'upper',
#' and \eqn{\frac{Pr(\mu_k=min\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k=min\{\mu_1,...,\mu_K\}){tp}}}} 
#' with \code{side} equals to 'lower', utilizing available data at each step.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_brar_unknown_var
#' @export sim_brar_unknown_var
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
#' @param au a vector of cut-off values in the final selection at the end of the trial,
#' with a length equal to the number of arms minus 1.
#' @param N2 maximal sample size for the trial.
#' @param tp tuning parameter. Some commonly used numbers are 0.5, 1 and n/2N.
#' @param armlabel a vector of treatment labels with an example of c(1, 2), where 1 and 2 describe
#' how each arm is labeled in a two-armed trial.
#' @param blocksize size of block used for equal randomization regarding participants in the 'initialization' period.
#' Recommend to be an even multiple of the number of total arms.
#' @param mean a vector of means in hypotheses, for example, as c(10,10) where 10 stands for the mean
#' in both groups. Another example is c(10,12) where 10 and 12 stand for the mean
#' for the control and the other treatment group, respectively.
#' @param sd a vector of standard deviations in hypotheses, for example, as c(2,2) where 2 stands for the standard deviation
#' in both groups. Another example is c(1,2) where 1 and 2 stand for the standard deviation
#' for the control and the other treatment group, respectively.
#' @param minstart a specified number of participants when one starts to check decision rules.
#' @param deltaa a vector of minimal effect expected to be observed for early futility stopping in
#' each arm is approximately \eqn{1\%}. The length of this parameter is \code{armn}-1.
#' @param tpp indicator of \code{tp} equals to n/2N. When \code{tp} is n/2N, \code{tpp} should be assigned 1. Default value is set to 0.
#' @param deltaa1 a vector of pre-specified minimal effect size expected to be observed at the final stage
#' for each arm. The length of this parameter is \code{armn}-1.
#' @param V01,a01,b01,m01 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 1 in the trial, which stands for the control.
#' @param V02,a02,b02,m02 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 2 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V03,a03,b03,m03 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 3 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V04,a04,b04,m04 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 4 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V05,a05,b05,m05 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 5 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return \code{sim_brar_unknown_var} returns an object of class "brar". An object of class "brar" is a list containing 
#' final decision, test statistics, the simulated data set and participants accrued for each arm 
#' at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' results. In the final decision, 'Superiorityfinal' refers to the selected arm, 
#' while 'Not Selected' indicates the arm stopped due to futility, and 'Control Selected' denotes the control arm chosen 
#' because other arms did not meet futility criteria before the final stage or were not deemed effective at the final stage. 
#' Note that before final stage of the trial, test statistics is calculated from \code{deltaa}, and test statistics is
#' calculated from \code{deltaa1} at the final stage.
#' @importFrom stats rnorm
#' @importFrom Rdpack reprompt
#' @examples
#' #sim_brar_unknown_var with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days under the null hypothesis,
#' #where mean=c(9.19/100,8.74/100,8.74/100), sd=c(0.009,0.009,0.009), tp=1 and 
#' #the minimal effect size is 0.
#' sim_brar_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length(vStartTime ),30,3)),enrollrate=0.1, N1=48,armn=3,au=c(0.85,0.85),
#' N2=480,tp=1,armlabel=c(1, 2,3),blocksize=6,mean=c(9.19/100,8.74/100,8.74/100),
#' sd=c(0.009,0.009,0.009), minstart=48,deltaa=c(-0.000325,-0.000325),
#' tpp=0,deltaa1=c(0,0),V01=1/2,a01=0.3,m01=9/100,b01=0.00001,side='lower')
#' @references 
#' \insertRef{Kevin2007}{RARtrials}


sim_brar_unknown_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,armn,au,N2,tp,armlabel,blocksize,
                               mean=0,sd=1,minstart,deltaa,tpp,deltaa1,V01,a01,b01,m01,
                               V02=V01,V03=V01,V04=V01,V05=V01,
                               a02=a01,a03=a01,a04=a01,a05=a01,
                               b02=b01,b03=b01,b04=b01,b05=b01,
                               m02=m01,m03=m01,m04=m01,m05=m01,side,...){

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
        data1[i,5]<-rnorm(1,mean=mean[j],
                          sd=sd[j])
      }
    }
  }

  armleft<-c(1:armn)
  decision<-rep(NA,armn )
  phi<-rep(NA,armn )
  stopp<-rep(NA,armn )
  V0=list(V01,V02,V03,V04,V05)
  a0=list(a01,a02,a03,a04,a05)
  b0=list(b01,b02,b03,b04,b05)
  m0=list(m01,m02,m03,m04,m05)
  
  for (jjj in minstart:N2){


    if (jjj>minstart){
      treat<-sample(armleft,size =1, prob = as.vector(pii))
      data1[jjj,4]<-treat
      data1[jjj,5]<-rnorm(1,mean=mean[treat],
                          sd=sd[treat])
    }

    if (jjj<N2){
      total<-sum (as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2]))
    }else if (jjj==N2){
      total<-N2
    }

    result<-vector("list",length(armleft))
    mat<-vector("list",armn)

    par<-vector("list",length(armleft))
    for (j in  1:length(armleft)){
      para<-list(V= V0[[armleft[j]]],a=a0[[armleft[j]]],
                 m=m0[[armleft[j]]],b=b0[[armleft[j]]])
     par[[j]]<-convert_gamma_to_chisq(para)
    }



    part<-vector("list",length(armleft))

    for (j in 1:length(armleft)) {

      if (total>0){
        if (jjj!=N2){
          data2<-matrix(data1[which(as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2])),],ncol=5)
        }else if (jjj==N2){
          data2<-data1
        }

        tot<-as.numeric(data2[which(data2[,4]==as.numeric(armleft[j])),5])

        if (identical(tot, numeric(0))){
          part[[armleft[j]]]<-par[[j]]
        }else{
          part[[armleft[j]]]<-update_par_nichisq(tot, par[[j]])
        }


      }else if (total==0){
        part[[armleft[j]]]<-par[[j]]
      }
    }

    if (length(armleft)>1){
      for (j in 1:length(armleft)){
        if (total>0){
          if (j>1){
            if (side=='lower'){
              result[[j]]<- pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa[armleft[j]-1],side='lower')

            }else if (side=='upper'){
              result[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa[armleft[j]-1],side='upper')

            }
          }else if (j==1){
            result[[1]]<-0
          }

        }else if (total==0 ){
          if (j>1){
            if (side=='lower'){
              result[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa[armleft[j]-1],side='lower')
            }else if (side=='upper'){
              result[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa[armleft[j]-1],side='upper')
            }

          }else if (j==1){
            result[[1]]<-0
          }
        }
      }
      aloo<-vector("list",length(armleft))
      aloo<-alofun_unk_var(mat=part,total=total,armleft=armleft,side=side)
    }


    if (jjj==N2){
      resultt<-vector("list",length(armleft))
      if (length(armleft)>1){
        for (j in 1:length(armleft)){
          if (total>0){
            if (j>1){
              if (side=='lower'){
                resultt[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='lower')
              }else if (side=='upper'){
                resultt[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='upper')
              }
            }else if (j==1){
              resultt[[1]]<-0
            }

          }else if (total==0 ){
            if (j>1){
              if (side=='lower'){
                resultt[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='lower')

              }else if (side=='upper'){
                resultt[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='upper')

              }

            }else if (j==1){
              resultt[[1]]<-0
            }
          }

        }
      }

    }


    posteriorp<-vector("list",length(armleft))
    for (j in 1:length(armleft)){
      if (total>0){
        if (j>1){
          if (side=='lower'){
            posteriorp[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='lower')

          }else if (side=='upper'){
            posteriorp[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='upper')

          }

        }else if (j==1){
          posteriorp[[j]]<-0
        }

      }else if (total==0 ){
        if (j>1){
          if (side=='lower'){
            posteriorp[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='lower')

          }else  if (side=='upper'){
            posteriorp[[j]]<-pgreater_NIX(part[[1]],part[[armleft[j]]],delta=deltaa1[armleft[j]-1],side='upper')

          }

        }else if (j==1){
          posteriorp[[j]]<-0
        }
      }
    }

    posteriorp1<-do.call(cbind,posteriorp)


    pii<-as.data.frame(do.call(cbind,aloo))
    colnames(pii)<-armleft
    if (jjj<N2){
      sim111<-do.call(cbind,result)
      colnames(sim111)<-armleft
    } else if (jjj==N2) {
      sim111t<-do.call(cbind,resultt)
      colnames(sim111t)<-armleft
    }

    if ( jjj<N2){


      for (k in 2:length(armleft)) {

        if (sim111[1,k]<0.01 & is.na(decision[armleft[k]])){
          decision[armleft[k]]<-'Futility'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]
        }
      }

      if ( 'Futility' %in% decision ){
        armleft<-armleft[! armleft %in%  which (decision %in%  'Futility')]
        pii<- pii[,sprintf("%s",armleft)]

      }

      if((length(armleft)==1 & length( which (decision %in%  'Futility'))==(armn-1))){
        stopp[ which (is.na(decision))]<-jjj
        data11<-data1[1:jjj,]
        if (is.na(decision[1])==TRUE) {
          decision[1]<-'Control Selected'
          phi[1]<-posteriorp1[1,1]
        }

        return(list(decision,phi,data11,stopp))
        break
      }


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

      for (k in 2:length(armleft)) {
        if (sim111t[1,k]>au[armleft[k]-1] & is.na(decision[armleft[k]])){
          decision[armleft[k]]<-'Superiorityfinal'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]

        }


        if (is.na(decision[armleft[k]] )) {

          decision[armleft[k]]<-'Not Selected'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]

        }
      }
      if (length( which (decision %in%  'Futility'))==(armn-1) &
          is.na(decision[1])){
        decision[1]<-'Control Selected'
        phi[1]<-posteriorp1[1,1]
      }
      nn<-rep(NA,armn)
      for (k in 1:armn) {
        nn[k]=nrow(data1[which(data1[,4]==k ),,drop=FALSE])
      }
      
      data11<-data1[1:N2,]
      output1<-list(decision[2:armn],phi[2:armn],data11,nn)
      class(output1)<-'brar'
      
     return(output1)
     
     # return(list(decision,phi,data11,nn))

    }
  }
}



#' @export 
print.brar<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}


