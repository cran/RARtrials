#' @title Calculate the Futility Stopping Probability for Continuous Endpoint with Unknown Variances Using a Normal-Inverse-Chi-Squared Distribution
#' @description Calculate the futility stopping probability in Bayesian response-adaptive randomization with
#' a control group using the Thall \eqn{\&} Wathen method for continuous outcomes with unknown variances. The prior distributions
#' follow Normal-Inverse-Chi-Squared (NIX) distributions and can be specified individually for each treatment group.
#' @details This function calculates the results of \eqn{Pr(\mu_k>\mu_{control}+\delta|data)} for \code{side} equals to
#' 'upper' and the results of \eqn{Pr(\mu_{control}>\mu_k+\delta|data)} for \code{side} equals to 'lower'.
#' The result indicates the posterior probability of stopping a treatment group due to futility around \eqn{1\%} in Bayesian
#' response-adaptive randomization with a control arm using Thall \eqn{\&} Wathen method, with accumulated results
#' during the conduct of trials. Parameters used in a Normal-Inverse-Gamma (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)})
#' distribution should be converted to parameters equivalent in a Normal-Inverse-Chi-Squared
#' (\eqn{(\mu,\sigma^2) \sim NIX(mean=\mu,effective sample size=\kappa,degrees of freedom=\nu,variance=\sigma^2/\kappa)})
#' distribution using \code{convert_gamma_to_chisq} before applying this function.
#' @aliases pgreater_NIX
#' @export pgreater_NIX
#' @param par1 current parameters including mu, kappa, nu, sigsq of a Normal-Inverse-Chi-Squared distribution from the control group.
#' @param par2 current parameters including mu, kappa, nu, sigsq of a Normal-Inverse-Chi-Squared distribution from the compared treatment group.
#' @param delta pre-specified minimal effect size expected to be observed between the control group and the compared treatment group.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to stats::integrate() (such as rel.tol) from this function.
#' @return a posterior probability of \eqn{Pr(\mu_k>\mu_{control}+\delta|data)} with \code{side} equals to 'upper';
#' a posterior probability of \eqn{Pr(\mu_{control}>\mu_k+\delta|data)} with \code{side} equals to 'lower'.
#' @importFrom stats pt
#' @importFrom stats dt
#' @examples
#' para<-list(V=1/2,a=0.5,m=9.1/100,b=0.00002)
#' par<-convert_gamma_to_chisq(para)
#' set.seed(123451)
#' y1<-rnorm(100,0.091,0.009)
#' par1<-update_par_nichisq(y1, par)
#' set.seed(123452)
#' y2<-rnorm(90,0.09,0.009)
#' par2<-update_par_nichisq(y2, par)
#' pgreater_NIX(par1=par1,par2=par2, side='upper')
#' pgreater_NIX(par1=par1,par2=par2, side='lower')
#' @references 
#' \insertRef{Wathen2017}{RARtrials}
#' \insertRef{Kevin2007}{RARtrials}

pgreater_NIX<-function(par1,par2,delta=0,side,...){
  f<-function(x) pnichisq_mu(x+delta,par2,side)*dnichisq_mu(x,par1)
  m<-par1$mu
  s<-sqrt(par1$sigsq/par1$kappa+par2$sigsq/par2$kappa)+abs(par1$mu-par2$mu)
  integrate(f,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
}

