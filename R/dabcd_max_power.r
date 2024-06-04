#' @title Allocation Probabilities Using Doubly Adaptive Biased Coin Design with Maximal Power Strategy for Binary Endpoint
#' @description \code{dabcd_max_power} can be used for doubly adaptive biased coin design with maximal power
#' strategy for binary outcomes, targeting generalized Neyman allocation and generalized RSIHR allocation. The return 
#' of this function is a vector of allocation probabilities to each arm, with the pre-specified number of participants in the trial.
#' @details The function simulates allocation probabilities for doubly adaptive biased coin design with maximal power strategy targeting
#' generalized Neyman allocation with 2-5 arms which is provided in \insertCite{Tymofyeyev2007}{RARtrials} or
#' generalized RSIHR allocation with 2-3 arms which is provided in \insertCite{Jeon2010}{RARtrials}, with modifications for typos
#' in \insertCite{Sabo2016}{RARtrials}. All of those methods are not smoothed. The output of this function is based on Hu \code{\&} Zhang's formula \insertCite{Hu2004}{RARtrials}.
#' With more than two armd the one-sided nominal level of each test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
#' @aliases dabcd_max_power
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export dabcd_max_power
#' @param NN a vector representing the number of participants with success results for each arm
#' estimated from the current data.
#' @param Ntotal1 a vector representing the total number of participants for each arm
#' estimated from the current data.
#' @param armn number of total arms in the trial.
#' @param BB the minimal allocation probability for each arm, which is within the
#' range of \eqn{[0,1/armn]}.
#' @param type allocation type, with choices from 'RSIHR' and 'Neyman'.
#' @param dabcd an indicator of whether to apply Hu & Zhang's formula (\insertCite{Hu2004}{RARtrials}), with choices from FALSE and TRUE.
#' TRUE represents allocation probabilities calculated using Hu & Zhang's formula;
#' FALSE represents allocation probabilities calculated before applying Hu & Zhang's formula.
#' Default value is set to FALSE.
#' @param gamma tuning parameter in Hu & Zhang's formula (\insertCite{Hu2004}{RARtrials}). When \code{dabcd}=FALSE, this parameter does not need
#' to be specified. Default value is set to 2.
#' @return A vector of allocation probabilities to each arm.
#' @examples
#' dabcd_max_power(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5,BB=0.2, type='Neyman')
#' dabcd_max_power(NN=c(54,67,85,63),Ntotal1=c(100,88,90,94),armn=4,BB=0.2, type='Neyman')
#' @references 
#' \insertRef{Hu2004}{RARtrials}
#' \insertRef{Tymofyeyev2007}{RARtrials}
#' \insertRef{Jeon2010}{RARtrials}
#' \insertRef{Sabo2016}{RARtrials}

dabcd_max_power<-function(NN,Ntotal1,armn,BB,type,dabcd=FALSE,gamma=2){

  NN<-NN+1
  Ntotal11<-Ntotal1+2
  p<-cbind(p=(NN/Ntotal11),arm=1:armn)
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
              rho[4]=sqrt(p5*q5)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
              rho[5]=sqrt(p5*q5)/(2*(sqrt(p1*q1)+sqrt(p5*q5)))
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
              rho[3]=sqrt(p5*q5)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
              rho[4]=sqrt(p5*q5)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
              rho[5]=sqrt(p5*q5)/(3*(sqrt(p1*q1)+sqrt(p5*q5)))
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


  if (dabcd==TRUE){
    alr<-rep(NA,armn)
    phi<-rep(NA,armn)
    for (k in 1:armn){
      phi[k]<-rho1[k]*((rho1[k]/(Ntotal1[k]/(sum(Ntotal1))))^gamma)
    }
    for (kk in 1:armn){
      alr[kk]<-phi[kk]/sum(phi)
    }
    return(alr)
  }else if (dabcd==FALSE){
    return(rho1)
  }

}

