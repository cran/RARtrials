#' @title Calculate the Futility Stopping Probability for Binary Endpoint with Beta Distribution
#' @description Calculate the futility stopping probability in Bayesian response-adaptive randomization with
#' a control group using the Thall \eqn{\&} Wathen method for binary outcomes. The conjugate prior distributions follow
#' Beta (\eqn{Beta(\alpha,\beta)}) distributions and can be specified individually for each treatment group.
#' @details This function calculates the results of \eqn{Pr(p_k>p_{control}+\delta|data)} for \code{side} equals to
#' 'upper' and the results of \eqn{Pr(p_{control}>p_k+\delta|data)} for \code{side} equals to 'lower'.
#' The result indicates the posterior probability of stopping a treatment group due to futility around \eqn{1\%} in Bayesian
#' response-adaptive randomization with a control arm using Thall \eqn{\&} Wathen method, with accumulated results
#' during the conduct of trials. 
#' @aliases pgreater_beta
#' @export pgreater_beta
#' @param a1,b1  \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)}, current estimated \eqn{\alpha} for the control group.
#' @param a2,b2  \eqn{\alpha} and \eqn{\beta} in \eqn{Beta(\alpha,\beta)}, current estimated \eqn{\alpha} for the treatment group which is compared to the control group.
#' @param delta  expected difference in success probabilities between the control group and the treatment group.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to stats::integrate() (such as rel.tol) from this function.
#' @return a posterior probability of \eqn{Pr(p_k>p_{control}+\delta|data)} with \code{side} equals to 'upper';
#' a posterior probability of \eqn{Pr(p_{control}>p_k+\delta|data)} with \code{side} equals to 'lower'.
#' @importFrom stats pbeta
#' @importFrom stats dbeta
#' @examples
#' pgreater_beta(a1=8, b1=10,a2=5, b2=19, delta=0.1, side='upper')
#' pgreater_beta(a1=65, b1=79,a2=58, b2=68, delta=0, side='lower')
#' @references 
#' \insertRef{Wathen2017}{RARtrials}


pgreater_beta<-function(a1,b1,a2,b2,delta,side,...){
  if (side=='upper'){
    f<-function(y){
      pbeta(y+delta,a2,b2,lower.tail=FALSE)*dbeta(y,a1,b1)
    }
  }else if (side=='lower'){
    f<-function(y){
      pbeta(y+delta,a2,b2,lower.tail=TRUE)*dbeta(y,a1,b1)
    }
  }
  integrate(f,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
}
