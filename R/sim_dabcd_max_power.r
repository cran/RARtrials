#' @title Simulate a Trial Using Doubly Adaptive Biased Coin Design with Maximal Power Strategy for Binary Endpoint
#' @description \code{sim_dabcd_max_power} can be used for doubly adaptive biased coin design with maximal power
#' strategy for binary outcomes, targeting generalized Neyman allocation and generalized RSIHR allocation. 
#' @details The function simulates a trial for doubly adaptive biased coin design with maximal power strategy targeting
#' generalized Neyman allocation with 2-5 arms which is provided in \insertCite{Tymofyeyev2007}{RARtrials} and
#' generalized RSIHR allocation with 2-3 arms which is provided in \insertCite{Jeon2010}{RARtrials}, with modifications for typos
#' in \insertCite{Sabo2016}{RARtrials}. All of those methods are not smoothed. The output of this function is based on Hu \code{\&} Zhang's formula \insertCite{Hu2004}{RARtrials}.
#' With more than two armd the one-sided nominal level of each test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
#' \code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
#' \code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
#' are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @aliases sim_dabcd_max_power
#' @export sim_dabcd_max_power
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
#' @param N1 number of participants with equal randomization in the burn-in period.
#' Recommend using 10 percent of the total sample size.
#' @param N2 maximal sample size for the trial.
#' @param armn number of total arms in the trial.
#' @param armlabel a vector of arm labels with an example of c(1, 2), where 1 and 2 describes
#' how each arm is labeled in a two-armed trial.
#' @param h a vector of success probabilities in hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the success probability
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probabilities
#' for the control and the treatment group, respectively.
#' @param BB the minimal allocation probabilities for each arm, which is within the
#' range of \eqn{[0,1/armn]}.
#' @param type allocation type, with choices from 'RSIHR' and 'Neyman'.
#' @param gamma tuning parameter in Hu & Zhang's formula. When dabcd=0, this parameter does not need
#' to be specified. Default value is set to 2.
#' @param alphaa the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @return \code{sim_dabcd_max_power} returns an object of class "dabcd". An object of class "dabcd" is a list containing 
#' final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected,
#' Z test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' result.
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @examples
#' sim_dabcd_max_power(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length( vStartTime ),30, 3)),enrollrate=0.9,N1=30,N2=300,armn=3,
#' armlabel=c(1,2,3),h=c(0.2,0.3,0.2),BB=0.1,type='Neyman',
#' side='upper')
#' sim_dabcd_max_power(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length( vStartTime ),60, 3)),enrollrate=0.1,N1=50,N2=500,armn=3,
#' armlabel=c(1,2,3),h=c(0.2,0.3,0.3),BB=0.15,type='RSIHR',
#' side='upper')
#'
#' @references 
#' \insertRef{Hu2004}{RARtrials}
#' \insertRef{Tymofyeyev2007}{RARtrials}
#' \insertRef{Jeon2010}{RARtrials}
#' \insertRef{Sabo2016}{RARtrials}

sim_dabcd_max_power<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,N2,
                              armn,armlabel,h,BB,type,gamma=2,alphaa=0.025,side){


  alr<-rep(NA_real_,armn)
  ncp<-rep(0,1)
  pr1<-rep(0,armn)

  popdat<-pop(Pats,nMax,enrollrate)
  data1<-startfun(popdat,TimeToOutcome=TimeToOutcome,blocksize=2*armn,N1=N1,armn=armn,armlabel=armlabel,N2=N2,h=h)
  low <- rep(0,armn)
  high <-rep(1,armn)

  for (i in N1:N2) {

    nt<-i

    if (i>N1){
      assigna<-sample(1:armn,size = 1,  prob = alr)
      data1[nt,4]<- assigna
      data1[nt,5]<-rbinom(1,size=1,prob=h[assigna])
    }
    if (nt<N2){
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[nt,2]))
    }else if (nt==N2){
      total1<-N2
    }

    if (total1>0){
      data2<-matrix(data1[which(as.numeric(data1[1:nt,3])<=as.numeric(data1[nt,2])),],ncol=5)
    }else if (total1==0){
      data2<-matrix(0,nrow=1,ncol=5)
    }

    dummy1<-matrix(NA,armn,2)
    for (m in 1:armn) {
      dummy1[m,1]<-sum(as.numeric(data2[data2[,4] %in% m,5]))
      dummy1[m,2]<-nrow(data2[data2[,4]== m,,drop=F])
    }

    NN<-dummy1[,1]+1
    Ntotal1<-dummy1[,2]+2


  p<-cbind(p=unname(unlist(NN/Ntotal1)),arm=1:armn)
  rank1<-cbind(p,rank=rank(-p[,'p'], ties.method = "min"),rankorder=rank(-p[,'p'],ties.method = "first"))
  rho<-rep(NA,armn)
  rho1<-rep(NA,armn)
  pk=rank1[rank1[,'rankorder'] == armn,'p']
  p1=rank1[rank1[,'rankorder'] == 1,'p']


  if (armn==2 & type=='Neyman'){
    rho1[2]<-sqrt(p[2,1]*(1-p[2,1]))/(sqrt(p[1,1]*(1-p[1,1]))+sqrt(p[2,1]*(1-p[2,1])))
    rho1[1]<-1-rho1[2]
  }
  if (armn==2 & type=='RSIHR'){
    rho1[2]<-sqrt(p[2,1])/sum(sqrt(p[1,1])+sqrt(p[2,1]))
    rho1[1]<-1-rho1[2]
  }

  if (armn==3 & type=='Neyman'){
    p1<-rank1[rank1[,'rankorder'] == 1,'p']
    p2<-rank1[rank1[,'rankorder'] == 2,'p']
    p3<-rank1[rank1[,'rankorder'] == 3,'p']
    q1<-1-p1
    q2<-1-p2
    q3<-1-p3
    if (length(unique(rank1[,'rank']))==3){
      ss=1
      gg=1

      por1<-NULL
      por2<-NULL
      recur<-c((ss+1):(armn-gg))

      for (j in 1:length(recur)) {
        por1[j]<- (pk*(1-pk))/((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
        por2[j]<-  ((rank1[rank1[,'rankorder'] == recur[j],'p'])-pk) /((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
      }

      alpha<-((sqrt(p1*(1-p1)))/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))*(sum(por1)-(armn-ss-gg))-
        ((sqrt(p1*(1-p1)*pk*(1-pk)))/(p1-pk))*sum(por2)

      B1tilt<-(1/(ss-alpha))*(sqrt(p1*(1-p1))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      BKtilt<-(1/(armn-ss+alpha))*(sqrt(pk*(1-pk))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      Btilt<-min(1/armn,B1tilt,BKtilt)


      if ((BB>Btilt) & (Btilt==B1tilt)){
        for (i in 1:(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-BB
        }
        for (ii in (armn-gg+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-(1-(armn-gg)*BB)/gg
        }

      }else if ((BB>Btilt) & (Btilt==BKtilt)) {

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1-(armn-ss)*BB)/ss
        }
        for (ii in (ss+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB
        }

      }else{

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1/ss)*(alpha*BB+sqrt(p1*(1-p1))/
                                                              (sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))
        }
        for ( ii in (ss+1):(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB
        }
        for ( iii in (armn-gg+1):(armn)) {
          rho[rank1[rank1[,'rankorder']==iii,'arm']]<-1/gg*(1-BB*(armn-ss-gg)-ss*rho[rank1[rank1[,'rankorder']==1,'arm']])
        }

      }
      rho1<-rho
    }else if (length(unique(rank1[,'rank']))==2){

      if (sum(rank1[,'rank']==3)==1){
        minn<-min(sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p3*q3))),sqrt(p3*q3)/(sqrt(p1*q1)+sqrt(p3*q3)),1/3 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p3*q3)))
          rho[2]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p3*q3)))
          rho[3]=sqrt(p3*q3)/(sqrt(p1*q1)+sqrt(p3*q3))
        }else if (BB > (sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p3*q3))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=1-2*BB
        }else if (BB > (sqrt(p3*q3)/(sqrt(p1*q1)+sqrt(p3*q3)))){
          rho[1]=(1-BB)/2
          rho[2]=(1-BB)/2
          rho[3]=BB
        }
      } else if (sum(rank1[,'rank']==2)==2){
        minn<-min(sqrt(p3*q3)/(2*(sqrt(p1*q1)+sqrt(p3*q3))),1/3 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(sqrt(p1*q1)+sqrt(p3*q3))
          rho[2]=sqrt(p3*q3)/(2*(sqrt(p1*q1)+sqrt(p3*q3)))
          rho[3]=sqrt(p3*q3)/(2*(sqrt(p1*q1)+sqrt(p3*q3)))
        }else if (BB > (sqrt(p3*q3)/(2*(sqrt(p1*q1)+sqrt(p3*q3))))){
          rho[1]=1-2*BB
          rho[2]=BB
          rho[3]=BB
        }

      }
      rho1[rank1[rank1[,'rankorder'] == 1,'arm']]<-   rho[1]
      rho1[rank1[rank1[,'rankorder'] == 2,'arm']]<-   rho[2]
      rho1[rank1[rank1[,'rankorder'] == 3,'arm']]<-   rho[3]
    }else if (length(unique(rank1[,'rank']))==1){
      rho1[1]<-1/3
      rho1[2]<-1/3
      rho1[3]<-1/3
    }

  }else if (armn==4 & type=='Neyman'){
    p1<-rank1[rank1[,'rankorder'] == 1,'p']
    p2<-rank1[rank1[,'rankorder'] == 2,'p']
    p3<-rank1[rank1[,'rankorder'] == 3,'p']
    p4<-rank1[rank1[,'rankorder'] == 4,'p']
    q1<-1-p1
    q2<-1-p2
    q3<-1-p3
    q4<-1-p4

    if (length(unique(rank1[,'rank'])) ==4|
        length(unique(rank1[,'rank']))==3){

      if (p1>p2 & (p2>p3 | p2==p3) & p3>p4 ){
        ss=1
        gg=1
      }else if (p1>p2 & p2>p3 & p3==p4){
        ss=1
        gg=2
      }else if (p1==p2 & p2>p3 & p3>p4){
        ss=2
        gg=1
      }

      por1<-NULL
      por2<-NULL
      recur<-c((ss+1):(armn-gg))

      for (j in 1:length(recur)) {
        por1[j]<- (pk*(1-pk))/((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
        por2[j]<-  ((rank1[rank1[,'rankorder'] == recur[j],'p'])-pk) /((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
      }

      alpha<-((sqrt(p1*(1-p1)))/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))*(sum(por1)-(armn-ss-gg))-
        ((sqrt(p1*(1-p1)*pk*(1-pk)))/(p1-pk))*sum(por2)

      B1tilt<-(1/(ss-alpha))*(sqrt(p1*(1-p1))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      BKtilt<-(1/(armn-ss+alpha))*(sqrt(pk*(1-pk))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      Btilt<-min(1/armn,B1tilt,BKtilt)


      if ((BB>Btilt) & (Btilt==B1tilt)){
        for (i in 1:(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-BB

        }
        for (ii in (armn-gg+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-(1-(armn-gg)*BB)/gg
        }

      }else if ((BB>Btilt) & (Btilt==BKtilt)) {

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1-(armn-ss)*BB)/ss
        }


        for (ii in (ss+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB
        }

      }else{

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1/ss)*(alpha*BB+sqrt(p1*(1-p1))/
                                                              (sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))
        }
        for ( ii in (ss+1):(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB

        }
        for ( iii in (armn-gg+1):(armn)) {
          rho[rank1[rank1[,'rankorder']==iii,'arm']]<-1/gg*(1-BB*(armn-ss-gg)-ss*rho[rank1[rank1[,'rankorder']==1,'arm']])
        }

      }

      rho1<-rho
    }else if (length(unique(rank1[,'rank']))==2){

      if (sum(rank1[,'rank']==1)==3){
        minn<-min(sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p4*q4))),sqrt(p4*q4)/(sqrt(p1*q1)+sqrt(p4*q4)),1/4 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[2]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[3]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[4]=sqrt(p4*q4)/(sqrt(p1*q1)+sqrt(p4*q4))
        }else if (BB > (sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p4*q4))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=BB
          rho[4]=1-3*BB
        }else if (BB > (sqrt(p4*q4)/(sqrt(p1*q1)+sqrt(p4*q4)))){
          rho[1]=(1-BB)/3
          rho[2]=(1-BB)/3
          rho[3]=(1-BB)/3
          rho[4]=BB
        }
      } else if (sum(rank1[,'rank']==2)==3){
        minn<-min(sqrt(p4*q4)/(3*(sqrt(p1*q1)+sqrt(p4*q4))),1/4 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(sqrt(p1*q1)+sqrt(p4*q4))
          rho[2]=sqrt(p4*q4)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[3]=sqrt(p4*q4)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[4]=sqrt(p4*q4)/(3*(sqrt(p1*q1)+sqrt(p4*q4)))
        }else if (BB > (sqrt(p4*q4)/(3*(sqrt(p1*q1)+sqrt(p4*q4))))){
          rho[1]=1-3*BB
          rho[2]=BB
          rho[3]=BB
          rho[4]=BB
        }

      }else if (sum(rank1[,'rank']==1)==2){
        minn<-min(sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p4*q4))),sqrt(p4*q4)/(2*(sqrt(p1*q1)+sqrt(p4*q4))),1/4 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[2]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[3]=sqrt(p4*q4)/(2*(sqrt(p1*q1)+sqrt(p4*q4)))
          rho[4]=sqrt(p4*q4)/(2*(sqrt(p1*q1)+sqrt(p4*q4)))
        }else if (BB > (sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p4*q4))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=(1-2*BB)/2
          rho[4]=(1-2*BB)/2
        }else if (BB > sqrt(p4*q4)/(2*(sqrt(p1*q1)+sqrt(p4*q4))) ){
          rho[1]=(1-2*BB)/2
          rho[2]=(1-2*BB)/2
          rho[3]=BB
          rho[4]=BB
        }

      }
      rho1[rank1[rank1[,'rankorder'] == 1,'arm']]<-   rho[1]
      rho1[rank1[rank1[,'rankorder'] == 2,'arm']]<-   rho[2]
      rho1[rank1[rank1[,'rankorder'] == 3,'arm']]<-   rho[3]
      rho1[rank1[rank1[,'rankorder'] == 4,'arm']]<-   rho[4]

    }else if (length(unique(rank1[,'rank']))==1){
      rho1[1]<-1/4
      rho1[2]<-1/4
      rho1[3]<-1/4
      rho1[4]<-1/4
    }



  }else if (armn==5 & type=='Neyman'){
    p1<-rank1[rank1[,'rankorder'] == 1,'p']
    p2<-rank1[rank1[,'rankorder'] == 2,'p']
    p3<-rank1[rank1[,'rankorder'] == 3,'p']
    p4<-rank1[rank1[,'rankorder'] == 4,'p']
    p5<-rank1[rank1[,'rankorder'] == 5,'p']
    q1<-1-p1
    q2<-1-p2
    q3<-1-p3
    q4<-1-p4
    q5<-1-p5
    if (length(unique(rank1[,'rank']))==5 |
        length(unique(rank1[,'rank']))==4 |
        length(unique(rank1[,'rank']))==3){


      if (p1>p2 & (p2>p3 |p2==p3) & (p3>p4 | p3==p4) & p4>p5){
        ss=1
        gg=1
      }else if (p1>p2 & (p2>p3 |p2==p3) & p3>p4 & p4==p5){
        ss=1
        gg=2
      }else if (p1==p2 & (p2>p3) & (p3>p4 | p3==p4) & p4>p5){
        ss=2
        gg=1
      }else if (p1==p2 & p2==p3 & p3>p4 & p4>p5){
        ss=3
        gg=1
      }else if (p1>p2 & p2>p3 & p3==p4 & p4==p5){
        ss=1
        gg=3
      }else if (p1==p2 & (p2>p3) & p3>p4 & p4==p5){
        ss=2
        gg=2
      }

      por1<-NULL
      por2<-NULL
      recur<-c((ss+1):(armn-gg))

      for (j in 1:length(recur)) {
        por1[j]<- (pk*(1-pk))/((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
        por2[j]<-  ((rank1[rank1[,'rankorder'] == recur[j],'p'])-pk) /((rank1[rank1[,'rankorder'] == recur[j],'p'])*(1-rank1[rank1[,'rankorder'] == recur[j],'p']))
      }

      alpha<-((sqrt(p1*(1-p1)))/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))*(sum(por1)-(armn-ss-gg))-
        ((sqrt(p1*(1-p1)*pk*(1-pk)))/(p1-pk))*sum(por2)

      B1tilt<-(1/(ss-alpha))*(sqrt(p1*(1-p1))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      BKtilt<-(1/(armn-ss+alpha))*(sqrt(pk*(1-pk))*(1/(sqrt(p1*(1-p1))+sqrt(pk*(1-pk)))))
      Btilt<-min(1/armn,B1tilt,BKtilt)


      if ((BB>Btilt) & (Btilt==B1tilt)){
        for (i in 1:(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-BB

        }
        for (ii in (armn-gg+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-(1-(armn-gg)*BB)/gg
        }

      }else if ((BB>Btilt) & (Btilt==BKtilt)) {

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1-(armn-ss)*BB)/ss
        }


        for (ii in (ss+1):armn) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB
        }

      }else{

        for (i in 1:ss) {
          rho[rank1[rank1[,'rankorder']==i,'arm']]<-(1/ss)*(alpha*BB+sqrt(p1*(1-p1))/
                                                              (sqrt(p1*(1-p1))+sqrt(pk*(1-pk))))
        }
        for ( ii in (ss+1):(armn-gg)) {
          rho[rank1[rank1[,'rankorder']==ii,'arm']]<-BB

        }
        for ( iii in (armn-gg+1):(armn)) {
          rho[rank1[rank1[,'rankorder']==iii,'arm']]<-1/gg*(1-BB*(armn-ss-gg)-ss*rho[rank1[rank1[,'rankorder']==1,'arm']])
        }
      }
      rho1<-rho
    }else if (length(unique(rank1[,'rank']))==2){

      if (sum(rank1[,'rank']==1)==4){
        minn<-min(sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5))),sqrt(p5*q5)/(sqrt(p1*q1)+sqrt(p5*q5)),1/5 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[2]=sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[3]=sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[4]=sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[5]=sqrt(p5*q5)/(sqrt(p1*q1)+sqrt(p5*q5))
        }else if (BB > (sqrt(p1*q1)/(4*(sqrt(p1*q1)+sqrt(p5*q5))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=BB
          rho[4]=BB
          rho[5]=1-4*BB
        }else if (BB > (sqrt(p5*q5)/(sqrt(p1*q1)+sqrt(p5*q5)))){
          rho[1]=(1-BB)/4
          rho[2]=(1-BB)/4
          rho[3]=(1-BB)/4
          rho[4]=(1-BB)/4
          rho[5]=BB
        }
      } else if (sum(rank1[,'rank']==2)==4){
        minn<-min(sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5))),1/5 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(sqrt(p1*q1)+sqrt(p5*q5))
          rho[2]=sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[3]=sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[4]=sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[5]=sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5)))
        }else if (BB > (sqrt(p5*q5)/(4*(sqrt(p1*q1)+sqrt(p5*q5))))){
          rho[1]=1-4*BB
          rho[2]=BB
          rho[3]=BB
          rho[4]=BB
          rho[5]=BB
        }

      }else if (sum(rank1[,'rank']==1)==3 ){
        minn<-min(sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5))),sqrt(p5*q5)/(2*(sqrt(p1*q1)+sqrt(p5*q5))),1/5 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[2]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[3]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[4]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[5]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
        }else if (BB > (sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=BB
          rho[4]=(1-3*BB)/2
          rho[5]=(1-3*BB)/2
        }else if (BB > sqrt(p5*q5)/(2*(sqrt(p1*q1)+sqrt(p5*q5))) ){
          rho[1]=(1-2*BB)/3
          rho[2]=(1-2*BB)/3
          rho[3]=(1-2*BB)/3
          rho[4]=BB
          rho[5]=BB
        }

      }else if (sum(rank1[,'rank']==1)==2 ){
        minn<-min(sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5))),sqrt(p5*q5)/(3*(sqrt(p1*q1)+sqrt(p5*q5))),1/5 )
        if (BB<=minn){
          rho[1]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[2]=sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[3]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[4]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
          rho[5]=sqrt(p1*q1)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
        }else if (BB > (sqrt(p1*q1)/(2*(sqrt(p1*q1)+sqrt(p5*q5))))){
          rho[1]=BB
          rho[2]=BB
          rho[3]=(1-2*BB)/3
          rho[4]=(1-2*BB)/3
          rho[5]=(1-2*BB)/3
        }else if (BB > sqrt(p5*q5)/(3*(sqrt(p1*q1)+sqrt(p5*q5))) ){
          rho[1]=(1-3*BB)/2
          rho[2]=(1-3*BB)/2
          rho[3]=BB
          rho[4]=BB
          rho[5]=BB
        }
      }
      rho1[rank1[rank1[,'rankorder'] == 1,'arm']]<-   rho[1]
      rho1[rank1[rank1[,'rankorder'] == 2,'arm']]<-   rho[2]
      rho1[rank1[rank1[,'rankorder'] == 3,'arm']]<-   rho[3]
      rho1[rank1[rank1[,'rankorder'] == 4,'arm']]<-   rho[4]
      rho1[rank1[rank1[,'rankorder'] == 5,'arm']]<-   rho[5]

    }else if (length(unique(rank1[,'rank']))==1){
      rho1[1]<-1/5
      rho1[2]<-1/5
      rho1[3]<-1/5
      rho1[4]<-1/5
      rho1[5]<-1/5
    }

  }

  if (armn==3 & type=='RSIHR'){
    p1=rank1[rank1[,'rankorder'] == 1,'p']
    p2=rank1[rank1[,'rankorder'] == 2,'p']
    p3=rank1[rank1[,'rankorder'] == 3,'p']
    q1=1-p1
    q2=1-p2
    q3=1-p3

    if (length(unique(rank1[,'rank']))==3){

      a= -(BB*q2-(BB-1)*q3) /(p1*q1)
      b= -(BB*(q3-q1))/(p2*q2)
      c= (BB*q2-(BB-1)*q1) /(p3*q3)
      d=sqrt(-a*b*(p1-p2)^2-a*c*(p1-p3)^2-b*c*(p2-p3)^2)
      l1=(a*(p1-p3)+b*(p2-p3)+d)/(p3*q3)
      l2=(b*(p1-p2)+c*(p1-p3)-d)/(p1*q1) +l1
      l3=(a*(p1-p2)-c*(p2-p3)+d)/(p2*q2) -l1

      rho[1]<-(l1+l3*BB)/l2
      rho[2]<-BB
      rho[3]<-1-BB-rho[1]

      if (rho[1]<=BB){
        rho[1]<-BB
        rho[2]<-BB
        rho[3]<-1-2*BB
      }else if (rho[3]<=BB){
        rho[1]<-1-2*BB
        rho[2]<-BB
        rho[3]<-BB
      }
    } else if (length(unique(rank1[,'rank']))==2 && sum(rank1[,'rank']==1)==2){

      if (BB<= min(sqrt(p1)/(2*(sqrt(p1)+sqrt(p3))),sqrt(p3)/(sqrt(p1)+sqrt(p3)), 1/3 )){
        rho[1]<-sqrt(p1)/(2*(sqrt(p1)+sqrt(p3)))
        rho[2]<-sqrt(p1)/(2*(sqrt(p1)+sqrt(p3)))
        rho[3]<-sqrt(p3)/(sqrt(p1)+sqrt(p3))
      } else if (BB> sqrt(p1)/(2*(sqrt(p1)+sqrt(p3)))){
        rho[1]<-BB
        rho[2]<-BB
        rho[3]<-1-2*BB
      } else if (BB>sqrt(p3)/(sqrt(p1)+sqrt(p3))){
        rho[1]<-(1-BB)/2
        rho[2]<-(1-BB)/2
        rho[3]<-BB
      }


    }else if (length(unique(rank1[,'rank']))==2 && sum(rank1[,'rank']==2)==2){

      if (BB<= min(sqrt(p3)/(2*(sqrt(p1)+sqrt(p3))),1/3 )){
        rho[1]<-sqrt(p1)/((sqrt(p1)+sqrt(p3)))
        rho[2]<-sqrt(p3)/(2*(sqrt(p1)+sqrt(p3)))
        rho[3]<-sqrt(p3)/(2*(sqrt(p1)+sqrt(p3)))
      } else if (BB> sqrt(p3)/(2*(sqrt(p1)+sqrt(p3)))){
        rho[1]<-1-2*BB
        rho[2]<-BB
        rho[3]<-BB
      }


    }else if (length(unique(rank1[,'rank']))==1){
      rho[1]<-1/3
      rho[2]<-1/3
      rho[3]<-1/3
    }


    rho1[rank1[rank1[,'rankorder'] == 1,'arm']]<-   rho[1]
    rho1[rank1[rank1[,'rankorder'] == 2,'arm']]<-   rho[2]
    rho1[rank1[rank1[,'rankorder'] == 3,'arm']]<-   rho[3]

  }

    alr<-rep(NA,armn)
    phi<-rep(NA,armn)
    if ( any(dummy1[,2]==0) ){
      phi<-rep(1/armn,armn)
    }else{
      for (k in 1:armn){
        phi[k]<-rho1[k]*((rho1[k]/(dummy1[k,2]/(sum(dummy1[,2]))))^gamma)
      }
    }

      for (kk in 1:armn){
        alr[kk]<-phi[kk]/sum(phi)
      }
    }

    pr<-vector("list",armn-1)
    phi<-vector("list",armn-1)
    dummy1<-matrix(NA,armn,2)
    for (m in 1:armn) {
      dummy1[m,1]<-sum(as.numeric(data1[data1[,4] %in% m,5]))
      dummy1[m,2]<-nrow(data1[data1[,4]==m,])
    }

    NN<-dummy1[,1]
    Ntotal1<-dummy1[,2]
    p<-cbind(p=unname(unlist(NN/Ntotal1)),arm=1:armn)

    for (l in 2:armn) {
      phi[[l-1]]<-(p[l,1]-p[1,1])/sqrt((p[l,1]*(1-p[l,1])/Ntotal1[l])+(p[1,1]*(1-p[1,1])/Ntotal1[1]))

      if (side=='upper'){
        if (phi[[l-1]]>=qnorm(1-alphaa/(armn-1))){
          pr[[l-1]]<-1 #success
        }else{
         pr[[l-1]]<-0
        }
      }else if (side=='lower'){
        if (phi[[l-1]]<=qnorm(alphaa/(armn-1))){
          pr[[l-1]]<-1 #success
        }else{
          pr[[l-1]]<-0
        }
      }

    }
    pr1<-do.call(cbind,pr)
    phi1<-do.call(cbind,phi)
   # return(list(pr1,phi1,data1))
    output1<-list(pr1,phi1,data1,Ntotal1)
    class(output1)<-'dabcd'
    
    
    return(output1)

}

#' @export 
print.dabcd<-function(x,...){
  cat("\nFinal Decision:\n",paste(x[[1]],sep=', ',collapse=', '),"\n")
  cat("\nTest Statistics:\n",paste(round(x[[2]],2),sep=', ',collapse=', '),"\n")
  cat("\nAccumulated Number of Participants in Each Arm:\n",paste(x[[4]],sep=', ',collapse=', '))
  invisible(x)
}