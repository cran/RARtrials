#' @title Posterior Probability that a Particular Arm is the Best for Binary Endpoint
#' @description Calculate posterior probability that a particular arm is the best in a trial using Bayesian response-adaptive randomization with
#' a control group (the Thall \eqn{\&} Wathen method). The conjugate prior distributions follow Beta (\eqn{Beta(\alpha,\beta)}) distributions
#' for binary outcomes in each arm and can be specified individually.
#' @details This function calculates the results of formula \eqn{Pr(p_k=max\{p_1,...,p_K\})} for
#' \code{side} equals to 'upper' and the results of formula \eqn{Pr(p_k=min\{p_1,...,p_K\})} for
#' \code{side} equals to 'lower'. This function returns the probability that the posterior probability of arm
#' \eqn{k} is maximal or minimal in trials with up to five arms.
#' @aliases pmax_beta
#' @export pmax_beta
#' @param armn number of arms in the trial with values up to 5. When \code{armn}=2,
#' only \code{a1} to \code{a2} and \code{b1} to \code{b2} need to be specified.
#' When \code{armn}=3, only \code{a1} to \code{a3} and  \code{b1} to \code{b3} need to be specified.
#' When \code{armn}=4, only \code{a1} to \code{a4} and  \code{b1} to \code{b4} need to be specified.
#' When \code{armn}=5, \code{a1} to \code{a5} and  \code{b1} to \code{b5} need to be specified.
#' @param a1,b1 \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)} for the arm to calculate the allocation probability of.
#' @param a2,b2 \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)} for one of the remaining arms.
#' @param a3,b3 \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)} for one of the remaining arms.
#' @param a4,b4 \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)} for one of the remaining arms.
#' @param a5,b5 \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)} for one of the remaining arms.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return a probability that a particular arm is the best in trials up to five arms.
#' @importFrom stats integrate
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats dnorm
#' @examples
#' pmax_beta(armn=5,a1=8,b1=10,a2=5,b2=19,a3=8,b3=21,
#' a4=6, b4=35, a5=15, b5=4, side='upper')
#' pmax_beta(armn=4,a1=56,b1=98,a2=25,b2=70,a3=87,b3=107,
#' a4=106, b4=202, side='lower')
#' pmax_beta(armn=3,a1=60,b1=46,a2=55,b2=46,a3=35,b3=36,side='upper')


pmax_beta<-function(armn,a1=NULL,b1=NULL,a2=NULL,b2=NULL,a3=NULL,b3=NULL,a4=NULL,b4=NULL,a5=NULL,b5=NULL,side,...){

  if (side=='upper'){

    if (armn==2){
      f1b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*dbeta(y,a1,b1)
      }
      integral<-integrate(f1b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==3){
      f2b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f2b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==4){
      f3b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          pbeta(y,a4,b4,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f3b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==5){
      f4b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          pbeta(y,a4,b4,lower.tail=TRUE)*pbeta(y,a5,b5,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f4b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }


  }else if (side=='lower'){

    if (armn==2){
      f5b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*dbeta(y,a1,b1)
      }
      integral<-integrate(f5b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==3){
      f6b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f6b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==4){
      f7b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          pbeta(y,a4,b4,lower.tail=FALSE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f7b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==5){
       f8b<-function(y){
          pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          pbeta(y,a4,b4,lower.tail=FALSE)*pbeta(y,a5,b5,lower.tail=FALSE)*
          dbeta(y,a1,b1)
       }
       integral<-integrate(f8b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }

  }
  


 return(integral)
}
