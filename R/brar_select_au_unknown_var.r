#' @title Select au in Bayesian Response-Adaptive Randomization with a Control Group for Continuous Endpoint with Unknown Variances
#' @description \code{brar_select_au_unknown_var} involves selecting au in Bayesian Response-Adaptive Randomization with a control group
#' for continuous endpoints with unknown variance in trials with two to five arms. The conjugate prior distributions follow
#' Normal-Inverse-Gamma (NIG) (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)})
#' distributions and can be specified individually for each arm.
#' @details This function generates a data set or a value in one iteration for selecting the appropriate au using Bayesian
#' response-adaptive randomization with a control group under null hypotheses with no delay and delayed scenarios.
#' The function can handle trials with up to 5 arms for continuous outcomes with unknown variances. This function uses the formula
#' \eqn{\frac{Pr(\mu_k=max\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k=max\{\mu_1,...,\mu_K\})^{tp}}}} with \code{side} equals to 'upper',
#' and \eqn{\frac{Pr(\mu_k=min\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k=min\{\mu_1,...,\mu_K\}){tp}}}} 
#' with \code{side} equals to 'lower', utilizing available data at each step.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases brar_select_au_unknown_var
#' @export brar_select_au_unknown_var
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
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param output control the output of brar_select_au_binary. If the user does not specify anything, the function returns
#' the entire dataset used to select the stopping boundary for each iteration. If the user specifies 'B', the function
#' only returns the selected stopping boundary for each iteration.
#' @param V01,a01,b01,m01 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 1 in the trial, which stands for the control.
#' @param V02,a02,b02,m02 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 2 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V03,a03,b03,m03 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 3 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V04,a04,b04,m04 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 4 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param V05,a05,b05,m05 prior parameters m, V, a, b in \eqn{NIG(V,m,a,b)} of arm 5 in the trial. Default value is set to \code{V01},
#' \code{a01}, \code{b01} and \code{m01}.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return A list of results generated from formula \eqn{Pr(\mu_k>\mu_{control}+\delta|data_{t-1})} at each step.
#' Note that before final stage of the trial, test statistics is calculated from \code{deltaa}, and test statistics is
#' calculated from \code{deltaa1} at the final stage.
#' @importFrom stats rnorm
#' @importFrom Rdpack reprompt
#' @examples
#' #brar_select_au_unknown_var with delayed responses follow a normal distribution with
#' #a mean of 60 days and a standard deviation of 3 days, where 
#' #mean=c((9.1/100+8.74/100+8.74/100)/3,(9.1/100+8.74/100+8.74/100)/3,
#' #(9.1/100+8.74/100+8.74/100)/3),sd=c(0.009,0.009,0.009),tp=1 and 
#' #the minimal effect size is 0. All arms have the same prior distributions.
#' set.seed(789)
#' stopbound1<-lapply(1:5,function(x){
#' brar_select_au_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length( vStartTime ),30, 3)), enrollrate=0.1, N1=48, armn=3, N2=480, tp=1,
#' armlabel=c(1,2,3), blocksize=6, mean=c((9.1/100+8.74/100+8.74/100)/3,
#' (9.1/100+8.74/100+8.74/100)/3,(9.1/100+8.74/100+8.74/100)/3) ,
#' sd=c(0.009,0.009,0.009),minstart=48, deltaa=c(-0.0003,-0.00035), tpp=0, 
#' deltaa1=c(0,0),V01=1/2,a01=0.3,m01=9/100,b01=0.00001,side='lower')
#'})
#'
#' simf<-list()
#' simf1<-list()
#' for (xx in 1:5){
#'  if (any(stopbound1[[xx]][48:479,2]<0.01)){
#'       simf[[xx]]<-NA
#'    }  else{
#'       simf[[xx]]<-stopbound1[[xx]][480,2]
#'  }
#'  if (any(stopbound1[[xx]][48:479,3]<0.01)){
#'       simf1[[xx]]<-NA
#'    }  else{
#'       simf1[[xx]]<-stopbound1[[xx]][480,3]
#'  }
#'}
#'simf2<-do.call(rbind,simf)
#'sum(is.na(simf2)) #1, achieve around 20% futility
#'simf3<-do.call(rbind,simf1)
#'sum(is.na(simf3)) #1, achieve around 20% futility
#'stopbound1a<-cbind(simf2,simf3)
#'stopbound1a[is.na(stopbound1a)] <- 0
#'sum(stopbound1a[,1]>0.85 | stopbound1a[,2]>0.85)/5 #0.2
#'#the selected stopping boundary is 0.85 with an overall lower one-sided type 
#'#I error of 0.2, based on 5 simulations. Because it is under the permutation null hypothesis,
#'#the selected deltaa should be an average of -0.0003 and -0.00035 which is -0.000325.
#'#It is recommended to conduct more simulations (i.e., 1000)  
#'#to obtain an accurate deltaa and au. As the simulation number increases, the
#'#choice of deltaa could be consistent for comparisons of each arm to the control.
#' @references 
#' \insertRef{Wathen2017}{RARtrials}

brar_select_au_unknown_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,armn,N2,tp,armlabel,blocksize,
                                     mean,sd, minstart,deltaa,tpp,deltaa1,side,output=NULL,V01,a01,b01,m01,
                                     V02=V01,V03=V01,V04=V01,V05=V01,
                                     a02=a01,a03=a01,a04=a01,a05=a01,
                                     b02=b01,b03=b01,b04=b01,b05=b01,
                                     m02=m01,m03=m01,m04=m01,m05=m01,...){

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
          data1[i,5]<-rnorm(1,mean=mean[j],sd=sd[j])
        }
      }
    }

    armleft<-c(1:armn)
    decision<-rep(NA,armn )
    phi<-rep(NA,armn )
    stopp<-rep(NA,armn )
    simout<-matrix(NA,nrow=N2,ncol=armn)
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


