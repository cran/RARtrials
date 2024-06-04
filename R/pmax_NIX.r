#' @title Posterior Probability that a Particular Arm is the Best for Continuous Endpoint with Unknown Variances 
#' @description Calculate posterior probability that a particular arm is the best in a trial using Bayesian response-adaptive randomization with
#' a control group (the Thall \eqn{\&} Wathen method). The conjugate prior distributions follow Normal-Inverse-Chi-Squared (NIX)
#' distributions for continuous outcomes with unknown variance in each arm and can be specified individually. 
#' @details This function calculates the results of formula \eqn{Pr(\mu_k=max\{\mu_1,...,\mu_k\})} for
#' \code{side} equals to 'upper' and the results of formula \eqn{Pr(\mu_k=min\{\mu_1,...,\mu_k\})} for
#' \code{side} equals to 'lower'. This function returns the probability that the posterior probability of arm
#' \eqn{k} is maximal or minimal in trials with up to five arms. Parameters used in a Normal-Inverse-Gamma 
#' (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)})
#' distribution should be converted to parameters equivalent in a Normal-Inverse-Chi-Squared
#' (\eqn{(\mu,\sigma^2) \sim NIX(mean=\mu,effective sample size=\kappa,degrees of freedom=\nu,variance=\sigma^2/\kappa)})
#' distribution using \code{convert_gamma_to_chisq} before applying this function.
#' @aliases pmax_NIX
#' @export pmax_NIX
#' @param armn number of arms in the trial with values up to 5. When \code{armn}=2,
#' only \code{par1} to \code{par2} need to be specified.
#' When \code{par1} to \code{par3} need to be specified.
#' When \code{par1} to \code{par4} need to be specified.
#' When \code{par1} to \code{par5} need to be specified.
#' @param par1 a vector of parameters including m, V, a, b for the arm with a Normal-Inverse-Chi-Squared prior
#' to calculate the allocation probability of.
#' @param par2 a vector of parameters including m, V, a, b for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par3 a vector of parameters including m, V, a, b for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par4 a vector of parameters including m, V, a, b for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par5 a vector of parameters including m, V, a, b for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return a probability that a particular arm is the best in trials up to five arms.
#' @importFrom stats integrate
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats dnorm
#' @examples
#' para<-list(V=1/2,a=0.8,m=9.1,b=1/2)
#' par<-convert_gamma_to_chisq(para)
#' set.seed(123451)
#' y1<-rnorm(100,9.1,1)
#' par11<-update_par_nichisq(y1, par)
#' set.seed(123452)
#' y2<-rnorm(90,9,1)
#' par22<-update_par_nichisq(y2, par)
#' set.seed(123453)
#' y3<-rnorm(110,8.92,1)
#' par33<-update_par_nichisq(y3, par)
#' y4<-rnorm(120,8.82,1)
#' par44<-update_par_nichisq(4, par)
#' pmax_NIX(armn=4,par1=par11,par2=par22,par3=par33,par4=par44,side='upper')
#' pmax_NIX(armn=4,par1=par11,par2=par22,par3=par33,par4=par44,side='lower')
#' 
#' para<-list(V=1/2,a=0.5,m=9.1/100,b=0.00002)
#' par<-convert_gamma_to_chisq(para)
#' set.seed(123451)
#' y1<-rnorm(100,0.091,0.009)
#' par11<-update_par_nichisq(y1, par)
#' set.seed(123452)
#' y2<-rnorm(90,0.09,0.009)
#' par22<-update_par_nichisq(y2, par)
#' set.seed(123453)
#' y3<-rnorm(110,0.0892,0.009)
#' par33<-update_par_nichisq(y3, par)
#' pmax_NIX(armn=3,par1=par11,par2=par22,par3=par33,side='upper')
#' pmax_NIX(armn=3,par1=par11,par2=par22,par3=par33,side='lower')


pmax_NIX<-function(armn,par1=NULL,par2=NULL,par3=NULL,par4=NULL,par5=NULL,side,...){
 
    if (armn==2){
      f1a1<-function(x) pnichisq_mu1(x,par2,side)*dnichisq_mu(x,par1)
      m<-par1$mu
      max_mu<-max(par1$mu,par2$mu)
      min_mu<-min(par1$mu,par2$mu)
      s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
      integral<-integrate(f1a1,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      
    }else if (armn==3){
      f1a2<-function(x) pnichisq_mu1(x,par2,side)*pnichisq_mu1(x,par3,side)*dnichisq_mu(x,par1)
      m<-par1$mu
      max_mu<-max(par1$mu,par2$mu,par3$mu)
      min_mu<-min(par1$mu,par2$mu,par3$mu)
      s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
      integral<-integrate(f1a2,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      
    }else if (armn==4){
      f1a3<-function(x) pnichisq_mu(x,par2,side)*pnichisq_mu(x,par3,side)*pnichisq_mu(x,par4,side)*dnichisq_mu(x,par1)
      m<-par1$mu
      max_mu<-max(par1$mu,par2$mu,par3$mu,par4$mu)
      min_mu<-min(par1$mu,par2$mu,par3$mu,par4$mu)
      s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
      integral<-integrate(f1a3,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      
    }else if (armn==5){
      f1a4<-function(x) pnichisq_mu(x,par2,side)*pnichisq_mu(x,par3,side)*pnichisq_mu(x,par4,side)*pnichisq_mu(x,par5,side)*dnichisq_mu(x,par1)
      m<-par1$mu
      max_mu<-max(par1$mu,par2$mu,par3$mu,par4$mu,par5$mu)
      min_mu<-min(par1$mu,par2$mu,par3$mu,par4$mu,par5$mu)
      s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
      integral<-integrate(f1a4,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
      
    }
  
  return(integral)
}
