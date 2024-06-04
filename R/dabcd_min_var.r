#' @title Allocation Probabilities Using Doubly Adaptive Biased Coin Design with Minimal Variance Strategy for Binary Endpoint
#' @description \code{dabcd_min_var} can be used for doubly adaptive biased coin design with minimal variance
#' strategy for binary outcomes, targeting generalized Neyman allocation and generalized RSIHR allocation. The return 
#' of this function is a vector of allocation probabilities to each arm, with the pre-specified number of participants in the trial.
#' @details The function simulates allocation probabilities for doubly adaptive biased coin design with minimal variance strategy targeting
#' generalized Neyman allocation and generalized RSIHR allocation with 2-5 arms. The output of this function is based on Hu \code{\&} Zhang's formula \insertCite{Hu2004}{RARtrials}.
#' With more than two armd the one-sided nominal level of each test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' @aliases dabcd_min_var
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export dabcd_min_var
#' @param NN a vector representing the number of participants with success results for each arm
#' estimated from the current data .
#' @param Ntotal1 a vector representing the total number of participants for each arm
#' estimated from the current data.
#' @param armn number of total arms in the trial.
#' @param type allocation type, with choices from 'RSIHR' and 'Neyman'.
#' @param dabcd an indicator of whether to apply Hu & Zhang's formula (\insertCite{Hu2004}{RARtrials}), with choices from FALSE and TRUE.
#' TRUE represents allocation probabilities calculated using Hu & Zhang's formula;
#' FALSE represents allocation probabilities calculated before applying Hu & Zhang's formula.
#' Default value is set to FALSE.
#' @param gamma tuning parameter in Hu & Zhang's formula (\insertCite{Hu2004}{RARtrials}). When \code{dabcd}=FALSE, this parameter does not need
#' to be specified. Default value is set to 2.
#' @return A vector of allocation probabilities to each arm.
#' @examples
#' dabcd_min_var(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5, type='Neyman')
#' dabcd_min_var(NN=c(54,67,85,63),Ntotal1=c(100,88,90,94),armn=4,type='RSIHR')
#' @references 
#' \insertRef{Hu2004}{RARtrials}

dabcd_min_var<-function(NN,Ntotal1,armn,type,dabcd=FALSE,gamma=2){

  NN<-NN+1
  Ntotal11<-Ntotal1+2
  p<-cbind(p=unname(unlist(NN/Ntotal11)),arm=1:armn)
  rho<-rep(NA,armn)
  rho1<-rep(NA,armn)

  if ( type=='Neyman'){
     for (k in 1:armn ){
       rho1[k]<-sqrt(p[k,1]*(1-p[k,1]))
     }
     for (k in 1:armn ){
       rho[k]<-rho1[k]/(sum(rho1))
     }
  }

  if (type=='RSIHR'){
     for (k in 1:armn ){
        rho1[k]<-sqrt(p[k,1])
     }
     for (k in 1:armn ){
        rho[k]<-rho1[k]/(sum(rho1))
     }
   }

  if (dabcd==FALSE){
     alr<-rep(NA,armn)
     phi<-rep(NA,armn)
     for (k in 1:armn){
       phi[k]<-rho[k]*((rho[k]/(Ntotal1[k]/(sum(Ntotal1)-1)))^gamma)
     }
     for (kk in 1:armn){
       alr[kk]<-phi[kk]/sum(phi)
     }
     return(alr)
  }else if (dabcd==TRUE){
     return(rho)
  }

}
