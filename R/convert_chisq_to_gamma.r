#' @title Convert parameters from a Normal-Inverse-Chi-Squared Distribution to a
#' Normal-Inverse-Gamma Distribution 
#' @description Convert parameters from a Normal-Inverse-Chi-Squared distribution to a
#' Normal-Inverse-Gamma distribution.
#' @details This function convert parameters from a Normal-Inverse-Chi-Squared
#' (\eqn{(\mu,\sigma^2) \sim NIX(mean=\mu,effective sample size=\kappa,degrees of freedom=\nu,variance=\sigma^2/\kappa)}) 
#' distribution to a Normal-Inverse-Gamma 
#' (\eqn{(\mu,\sigma^2) \sim NIG(mean=m,variance=V \times \sigma^2,shape=a,rate=b)}) 
#' distribution.
#' @aliases convert_chisq_to_gamma
#' @export convert_chisq_to_gamma
#' @param cpar a list of parameters including mu, kappa, nu, sigsq from a Normal-Chi-Squared distribution. 
#' @return a list of parameters including m, V, a, b from a Normal-Inverse-Gamma distribution.
#' @examples
#' convert_chisq_to_gamma(list(mu=0.091,kappa=2,nu=1,sigsq=4e-05))
#' @references 
#' \insertRef{Kevin2007}{RARtrials}

convert_chisq_to_gamma<-function(cpar){
  list(
    m=cpar$mu,
    V=1/cpar$kappa,
    a=cpar$nu/2,
    b=cpar$nu*cpar$sigsq/2
  )
}



