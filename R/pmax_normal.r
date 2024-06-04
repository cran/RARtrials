#' @title Posterior Probability that a Particular Arm is the Best for Continuous Endpoint with Known Variances
#' @description Calculate posterior probability that a particular arm is the best in a trial using Bayesian response-adaptive randomization with
#' a control group (the Thall \eqn{\&} Wathen method). The conjugate prior distributions follow Normal (\eqn{N(mean,sd)}) distributions for 
#' continuous outcomes with known variance in each arm and can be specified individually.
#' @details This function calculates the results of formula \eqn{Pr(\mu_k=max\{\mu_1,...,\mu_K\})} for
#' \code{side} equals to 'upper' and the results of formula \eqn{Pr(\mu_k=min\{\mu_1,...,\mu_K\})} for
#' \code{side} equals to 'lower'. This function returns the probability that the posterior probability of arm
#' \eqn{k} is maximal or minimal in trials with up to five arms.
#' @aliases pmax_normal
#' @export pmax_normal
#' @param armn number of arms in the trial with values up to 5. When \code{armn}=2,
#' only \code{mean1} to \code{mean2} and \code{sd1} to \code{sd2} need to be specified.
#' When \code{armn}=3, only \code{mean1} to \code{mean3} and \code{sd1} to \code{sd3} need to be specified.
#' When \code{armn}=4, only \code{mean1} to \code{mean4} and \code{sd1} to \code{sd4} need to be specified.
#' When \code{armn}=5, \code{mean1} to \code{mean5} and \code{sd1} to \code{sd5} need to be specified.
#' @param mean1,sd1 mean and sd in Normal(mean,sd) for the arm to calculate the allocation probability of.
#' @param mean2,sd2 mean and sd in Normal(mean,sd) for one of the remaining arms.
#' @param mean3,sd3 mean and sd in Normal(mean,sd) for one of the remaining arms.
#' @param mean4,sd4 mean and sd in Normal(mean,sd) for one of the remaining arms.
#' @param mean5,sd5 mean and sd in Normal(mean,sd) for one of the remaining arms.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return a probability that a particular arm is the best in trials up to five arms.
#' @importFrom stats integrate
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats dnorm
#' @examples
#' pmax_normal(armn=5,mean1=0.8,sd1=0.2,mean2=0.5,sd2=0.1,mean3=0.8,
#' sd3=0.5,mean4=0.6,sd4=0.2,mean5=0.6,sd5=0.2,side='upper')
#' pmax_normal(armn=4,mean1=8,sd1=2,mean2=8.5,sd2=2,mean3=8.3,
#' sd3=1.8,mean4=8.7,sd4=2,side='lower')
#' pmax_normal(armn=3,mean1=80,sd1=20,mean2=50,sd2=10,mean3=80,
#' sd3=15,side='upper')

pmax_normal<-function(armn,mean1=NULL,sd1=NULL,mean2=NULL,sd2=NULL,mean3=NULL,sd3=NULL,mean4=NULL,
                      sd4=NULL,mean5=NULL,sd5=NULL,side,...){

    
    if (side=='upper'){
      
      if (armn==2){
        f1c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f1c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==3){
        f2c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f2c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==4){
        f3c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
            pnorm(y,mean=mean4,sd=sd4,lower.tail=T)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f3c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==5){
        f4c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
            pnorm(y,mean=mean4,sd=sd4,lower.tail=T)*
            pnorm(y,mean=mean5,sd=sd5,lower.tail=T)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f4c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      }
      
    }else if (side=='lower'){
      
      if (armn==2){
        f5c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f5c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==3){
        f6c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f6c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==4){
        f7c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
            pnorm(y,mean=mean4,sd=sd4,lower.tail=F)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f7c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      } else if (armn==5){
        f8c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
            pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
            pnorm(y,mean=mean4,sd=sd4,lower.tail=F)*
            pnorm(y,mean=mean5,sd=sd5,lower.tail=F)*
            dnorm(y,mean=mean1,sd=sd1)
        }
        integral<-integrate(f8c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      }
      
    }
    
  return(integral)
}
