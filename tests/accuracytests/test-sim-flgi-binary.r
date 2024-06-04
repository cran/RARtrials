
########The set up of the code is to match the Table 1 & Table 3 results from Villar, S. S., Wason, J., & Bowden, J. (2015).
######## Response-adaptive randomization for multi-arm clinical trials using the forward looking Gittins index rule.
######## Biometrics, 71(4), 969–978. https://doi.org/10.1111/biom.12337


options(scipen=999)
SimulateAMonthOfAccrualTimes <- function( dPatsPerMonth , dStartMonth )
{
  nQtyPats    <- 1.2 *qpois(0.9999,dPatsPerMonth)
  vTimes      <- cumsum( rexp( nQtyPats, dPatsPerMonth ) )
  vTimes      <- vTimes[ vTimes < 1 ]
  vTimes      <- vTimes + dStartMonth
  return( vTimes )
}


SimulateArrivalTimes <- function( vPatsPerMonth, nMaxQtyPats )
{
  vTimes <- c()
  if( length( vPatsPerMonth ) == 1 )
  {
    vTimes <- cumsum(rexp(nMaxQtyPats ,vPatsPerMonth))
  }
  else
  {
    dStartMonth <- 0
    nMonth     <- 1
    while( length( vTimes ) < nMaxQtyPats  )
    {

      vTimes      <- c( vTimes, SimulateAMonthOfAccrualTimes( vPatsPerMonth[ nMonth ], dStartMonth ))
      dStartMonth <- dStartMonth + 1

      if( nMonth < length( vPatsPerMonth ) )
        nMonth <- nMonth +  1
    }
    vTimes <- vTimes[ 1:nMaxQtyPats ]
  }
  return( vTimes )
}

SimulateOutcomeObservedTime <- function( vStartTime )
{
  vTimeToOutcome <- 0
  vObsTime <- vStartTime  + vTimeToOutcome
  return( vObsTime )
}


allocation_probabilities<-function(tt,data1,K1,I0,block,noRuns,rule){

  index<-matrix(0,nrow=K1,1)
  selected<-matrix(0,nrow=noRuns,block)
  prob<-matrix(0,K1,block)

  for (j in 1:noRuns) {
    n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
    s=matrix(I0[,1]+1,nrow=nrow(I0),1)
    f=matrix(I0[,2]+1,nrow=nrow(I0),1)

    for (t1 in 0: (block-1)){
      for (k in 1:K1){
        index[k,1]=Gittins[f[k,1],s[k,1]]
      }
      max_index=max(index)
      kmax=min(match(max(index), index))
      idx= which(as.vector(index)==max(index))
      if (length(idx)>1){
        posi=sample(1:length(idx),1)
        kmax=idx[posi]
        max_index=as.vector(index)[kmax]
      }
      selected[j,t1+1]=kmax
      snext=s
      fnext=f
      nnext=n
      if (rule=='FLGI PM' | rule=='Controlled FLGI' ){#| rule=='Controlled FLGI modified'
        probSuc_kmax=s[kmax,1]/(s[kmax,1]+f[kmax,1])
        if(runif(1)<=probSuc_kmax){
          Pos=1
        }else{
          Pos=0
        }
      }else if (rule=='FLGI PD'){
        probSuc_kmax=rbeta(1,s[kmax,1],f[kmax,1])
        Pos= rbinom(1,1,probSuc_kmax)
      }
      nnext[kmax,1]=n[kmax,1]+1

      data1[tt*block+t1+1,4]=kmax
      data1[tt*block+t1+1,5]=Pos
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2]))

      for (k in 1:K1){
        if (total1>0){

          dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2])),],ncol=5)
          snext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==1,,drop=F])+2
          fnext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==0,,drop=F])+2

        }else if (total1==0){
          snext[k,1]=s[k,1]
          fnext[k,1]=f[k,1]
        }
      }

      s=snext
      f=fnext
      n=nnext
    }
  }

  for (i in 1:block){
    for(k in 1:K1){
      prob[k,i]=sum(selected[,i]==k)/noRuns
    }
  }

  allocation_probabilities=rowMeans(prob)

  return(allocation_probabilities)
}

allocation_probabilities1<-function(tt,data1,K1,I0,block,noRuns,rule){

  index<-matrix(0,nrow=K1,1)
  selected<-matrix(0,nrow=noRuns,block)
  prob<-matrix(0,K1,block)

  for (j in 1:noRuns) {
    n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
    s=matrix(I0[,1]+1,nrow=nrow(I0),1)
    f=matrix(I0[,2]+1,nrow=nrow(I0),1)

    for (t1 in 0: (block-1)){
      for (k in 1:K1){
        index[k,1]=Gittins[f[k,1],s[k,1]]
      }
      max_index=max(index)
      kmax=min(match(max(index), index))
      idx= which(as.vector(index)==max(index))
      if (length(idx)>1){
        posi=sample(1:length(idx),1)
        kmax=idx[posi]
        max_index=as.vector(index)[kmax]
      }
      selected[j,t1+1]=kmax
      snext=s
      fnext=f
      nnext=n
      if (rule=='FLGI PM'|rule=='Controlled FLGI'){
        probSuc_kmax=s[kmax,1]/(s[kmax,1]+f[kmax,1])
        if(runif(1)<=probSuc_kmax){
          Pos=1
        }else{
          Pos=0
        }
      }else if (rule=='FLGI PD'){
        probSuc_kmax=rbeta(1,s[kmax,1],f[kmax,1])
        Pos= rbinom(1,1,probSuc_kmax)
      }
      nnext[kmax,1]=n[kmax,1]+1
      data1[tt*block+t1+1,4]=kmax+1
      data1[tt*block+t1+1,5]=Pos
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2]))

      for (k in 1:K1){
        if (total1>0){
          dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[tt*block+t1+1,2])),],ncol=5)
          snext[k,1]=nrow(dataa[dataa[,4]==k+1 & dataa[,5]==1,,drop=F])+2
          fnext[k,1]=nrow(dataa[dataa[,4]==k+1 & dataa[,5]==0,,drop=F])+2

        }else if (total1==0){
          snext[k,1]=s[k,1]
          fnext[k,1]=f[k,1]
        }
      }
      s=snext
      f=fnext
      n=nnext
    }
  }

  for (i in 1:block){
    for(k in 1:K1){
      prob[k,i]=sum(selected[,i]==k)/noRuns
    }
  }

  allocation_probabilities=rowMeans(prob)

  return(allocation_probabilities)
}



simula_ocs_adapt_random<-function(jj,I0,K,crit,noRuns2,noRuns,Tsize,ptrue,block,criteria,rule,options){
  index<-matrix(0,nrow=K,1)
  thisSum<-matrix(0,nrow=noRuns,1)
  phat<-matrix(0,nrow=noRuns,K)
  sigmahat<-matrix(0,nrow=noRuns,K)
  # y<-matrix(0,nrow=noRuns,Tsize)
  ns<-matrix(0,nrow=noRuns,K)
  sn<-matrix(0,nrow=noRuns,K)
  zs1<-matrix(0,nrow=noRuns,K-1)
  pv<-matrix(0,nrow=noRuns,K-1)
  pv0<-matrix(0,nrow=noRuns,K-1)
  zp1<-matrix(0,nrow=noRuns,K-1)
  pvalues1<-matrix(0,nrow=noRuns,K-1)
  pvaluep1<-matrix(0,nrow=noRuns,K-1)
  pbetter<-matrix(0,nrow=noRuns,K-1)
  selected<-matrix(0,nrow=Tsize,noRuns)
  ntreat<-matrix(0,nrow=noRuns,K)
  ap<-matrix(0,nrow=noRuns,K-1)
  nbetter<-matrix(0,nrow=noRuns,K)
  besteffect<-max(ptrue)
  besttreat<-min(which(as.vector(ptrue) %in% max(ptrue)))


  for (j in 1:noRuns){
    vStartTime<-sort(sim11[[jj]][[j]][[3]][1:Tsize], decreasing = FALSE)
    vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime)
    data1<-matrix(NA,nrow=Tsize,ncol=5)
    data1[,1]<-1:Tsize
    data1[,2]<-vStartTime
    data1[,3]<-vOutcomeTime

    n=matrix(rowSums(I0)+2,nrow=nrow(I0),1)
    s=matrix(I0[,1]+1,nrow=nrow(I0),1)
    f=matrix(I0[,2]+1,nrow=nrow(I0),1)

    for (t in 0:((Tsize/block)-1)){
      alp=allocation_probabilities(tt=t,data1=data1,I0=cbind(s-1,f-1),block=block,noRuns=noRuns2,K1=K,rule=rule)

      if (rule=='FLGI PM' |rule=='FLGI PD'  ){
        crit=0.056
        if (K>2 & block>1){
          crit=0.05
        }else if (K>2 & block==1){
          crit=0.06
        }
        crit=(K-1)*crit
      }
      if (rule=='Controlled FLGI'  ){
        alp[1]=1/(K-1)
        elp_e=allocation_probabilities1(tt=t,data1=data1,I0=cbind(s[2:K,]-1,f[2:K,]-1),block=block,noRuns=noRuns2,K1=K-1,rule='FLGI PM')
        c=alp[1]+sum(elp_e)
        alp=(1/c)*c(alp[1],elp_e)
      }


      alp=cumsum(c(0,alp))

      snext=s
      fnext=f
      nnext=n

      Pob<-rep(0,block)
      Pos<-rep(0,block)
      for (p in 1:block){
        Pob[p]<-runif(1)
        for (k in 1:K){
          if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
            nnext[k]=n[k]+1
            if (runif(1)<=ptrue[k]){
              Pos[p]=1
            }else{
              Pos[p]=0
            }
            data1[t*block+p,4]=k
            data1[t*block+p,5]=Pos[p]
          }
        }
        total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[t*block+p,2]))

        for (k in 1:K){
          if (total1>0){
            dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[t*block+p,2])),],ncol=5)
            snext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==1,,drop=F])+2
            fnext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==0,,drop=F])+2
          }else if (total1==0){
            snext[k,1]=s[k,1]
            fnext[k,1]=f[k,1]
          }
        
       }

        s=snext
        f=fnext
        n=nnext
      }
      thisSum[j,1]=thisSum[j]+sum(Pos)

    }

    if (floor(Tsize/block)*block!=Tsize){
      Pob<-rep(0,Tsize-floor(Tsize/block)*block)
      Posi<-rep(0,Tsize-floor(Tsize/block)*block)
      for (p in 1:(Tsize-floor(Tsize/block)*block)){
        Pob[p]<-runif(1)
        for (k in 1:K){
          if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
            nnext[k]=n[k]+1
            if (runif(1)<=ptrue[k]){
              Posi[p]=1
            }else {
              Posi[p]=0
            }
            data1[floor(Tsize/block)*block+p,4]=k
            data1[floor(Tsize/block)*block+p,5]=Posi[p]
          }
        }
        total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[floor(Tsize/block)*block+p,2]))
        
        for (k in 1:K){
          if (total1>0){
            dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[floor(Tsize/block)*block+p,2])),],ncol=5)
            snext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==1,,drop=F])+2
            fnext[k,1]=nrow(dataa[dataa[,4]==k & dataa[,5]==0,,drop=F])+2
          }else if (total1==0){
            snext[k,1]=s[k,1]
            fnext[k,1]=f[k,1]
          }
        }

        s=snext
        f=fnext
        n=nnext
      }
      thisSum[j]=thisSum[j]+sum(Posi)

    }


    for (k in 1:K){
      s[k,1]=nrow(data1[data1[,4]==k & data1[,5]==1,,drop=F])+2
      f[k,1]=nrow(data1[data1[,4]==k & data1[,5]==0,,drop=F])+2
      n[k,1]=nrow(data1[data1[,4]==k ,,drop=F])+4
    }
    ns[j,]=n-2
    sn[j,]=s-1
    phat[j,]=(s-1)/(n-2)
    sigmahat[j,]=(phat[j,]*(1-phat[j,]))/ns[j,]

    sigma<-matrix(0,K-1,K-1)
    sigmat<-matrix(0,K-1,K-1)
    pc<-matrix(0,j,K-1)

    for (k in 1:(K-1)){

      pv[j,k]=phyper(sn[j,k+1]-1,sn[j,1]+sn[j,k+1],ns[j,1]+ns[j,k+1]-sn[j,1]-sn[j,k+1],ns[j,k+1])
      pv0[j,k]=phyper(sn[j,1]-1,sn[j,1]+sn[j,k+1],ns[j,k+1]+ns[j,1]-sn[j,1]-sn[j,k+1],ns[j,1])
      zs1[j,k]=(phat[j,k+1]-phat[j,1])/sqrt(sigmahat[j,1]+sigmahat[j,k+1])

    }

    ntreat[j,]=ns[j,]/sum(ns[j,])
    nbetter[j,]=n[besttreat,]/sum(n)
  }


  if (criteria==1){
    ap=(matrix(1,noRuns,K-1)-pv)-(crit/(K-1))*(matrix(1,noRuns,K-1))
    ap0=(matrix(1,noRuns,K-1)-pv0)-(crit/(K-1))*(matrix(1,noRuns,K-1))
    crit0=crit/2
    ap_1=(matrix(1,noRuns,K-1)-pv)-(crit0/(K-1))*(matrix(1,noRuns,K-1))
    ap_0=pv-(crit0/(K-1))*(matrix(1,noRuns,K-1))

  }else if (criteria==0){
    crit0=0.05/2
    ap_0=(matrix(1,noRuns,K-1)-pnorm(zs1))-(crit0/(K-1))*(matrix(1,noRuns,K-1))
    ap_1=pnorm(zs1)-(crit0/(K-1))*(matrix(1,noRuns,K-1))
    ap= (matrix(1,noRuns,K-1) -pnorm(zs1))-(crit/(K-1))*(matrix(1,noRuns,K-1))
    if (max(ptrue)!=ptrue[1]){
      ap0=(matrix(1,noRuns,K-1)-pnorm(zs1))-(crit/(K-1))*(matrix(1,noRuns,K-1))
      ap1=pnorm(zs1)-(crit/(K-1))*(matrix(1,noRuns,K-1))

    }else if (max(ptrue)==ptrue[1] & sum(ptrue==max(ptrue))<K){
      ap0=(matrix(1,noRuns,K-1)-pnorm(zs1))-(crit/(K-1))*(matrix(1,noRuns,K-1))
      ap1=pnorm(zs1)-(crit/(K-1))*(matrix(1,noRuns,K-1))
    }

  }

  b<-matrix(0,nrow=noRuns,1)
  b1<-matrix(0,nrow=noRuns,1)
  b2<-matrix(0,nrow=noRuns,1)
  p<-matrix(0,nrow=noRuns,1)
  e<-matrix(0,nrow=noRuns,1)
  for (h in (1:noRuns)){
    if (criteria==0){#something we do not understand in orginal code and has been updated
      b1[h,]= ifelse(sum(ap[h,]<0)>0,1,0)
      b[h,]= ifelse(sum(ap_0[h,]<0)>0 |sum(ap_1[h,]<0)>0 ,1,0)
      p[h,]<-ifelse(sum(b[h,])==0,1,0)
      e[h,]<-1-p[h,]
    }else{
      b1[h,]=ifelse(sum(ap[h,]<0)>0,1,0)
      b[h,]=ifelse(sum(ap_0[h,]<0)>0 |sum(ap_1[h,]<0)>0 ,1,0)
      p[h,]<-ifelse(sum(b[h,])==0,1,0)
      e[h,]<-1-p[h,]
    }
  }

  globalerror= mean(e)
  meanoc= globalerror
  meanoc_H= mean(b1)
  meanw= 'na'
  pw=matrix(0,nrow=noRuns,1)

  if (K>2 & sum(abs(ptrue-ptrue[1]*rep(1,K)))>0){
    ini=abs(ptrue-ptrue[1]*rep(1,K))
    sub=which(ini %in% max(ini))
    if (criteria==0 & max(ptrue)==ptrue[1] & sum(ptrue==max(ptrue))<K){
      subw=min(which(ptrue %in% max(ptrue[2:K]))-1)
      meanoc=mean(ap1[,subw]<0)
      subw=min(which(ptrue %in% min(ptrue[2:K]))-1)
      meanw=mean(ap1[,sub-1]<0)
      meanoc_H=mean(b1)
      #mean(ap1<0)
    }else if (criteria==0 & max(ptrue)!=ptrue[1]){
      meanoc=mean(ap[,sub-1]<0)
      meanoc_H=mean(b1)
      meanw= 'na'
    }

    if (criteria!=0 & max(ptrue)!=ptrue[1]){
      meanoc=mean(ap[,sub-1]<0)
      mean_oc=mean(ap<0)
      meanoc_H=mean(b1)
      meanw= 'na'

    }else if (criteria!=0 & max(ptrue)==ptrue[1] ){ #& sum(ptrue==max(ptrue))< K
      subw=min(which(ptrue %in% max(ptrue[2:K]))-1)
      meanoc=mean(ap0[,subw]<0)
      meanoc_H=mean(sum(ap0<0)>0|sum(ap<0)>0)
      meanoc_H=mean(b1)
      subw=min(which(ptrue %in% min(ptrue[2:K]))-1)
      meanw=mean(ap0[,sub-1]<0)
    }
  }

  decision<-matrix(0,1,K)
  for (k in 2:K+1){
    decision[1,k-1]= sum(pw==(k-2))/noRuns;
  }
  if  (K>4){
    decision='na'
  }

  #### If options=1, then selection is done based on GI criteria
  indexa<-matrix(0,noRuns,K)
  if (options==1){
    for (k in 1:K){
      if (rule %in% c('FLGI PM', 'Controlled Gittins','FLGI PD')){
        for (h in 1:noRuns){
          indexa[h,k] = Gittins[ns[h,k]-sn[h,k]+2,sn[h,k]+1]#Gittins[sn[h,k]+1,ns[h,k]-sn[h,k]+2]
        }
        pw=max.col(indexa) #need calculate the column maximal
        pw=pw-1
      }
    }
    for (k in 2:K+1){
      decision[1,k-1]= sum(pw==(k-2))/noRuns;
    }
  }

  meanpbest=mean(nbetter)
  av_mpb=sd(nbetter)
  av_value = mean(thisSum)
  av_std = sd(thisSum)
  mnum= mean(ns)

  return(list(av_value,av_std,meanoc,meanpbest,av_mpb,meanw,meanoc_H))
}


####################Run before other code ############
pop<-function(vPatsPerMonth,nMaxQtyPats,enrollrate1){
  populationtotal<-SimulateArrivalTimes (vPatsPerMonth, nMaxQtyPats)

  vStartTime1<-rbinom(nMaxQtyPats,size=1,enrollrate1)
  vStartTime2<- cbind(vStartTime1,populationtotal)
  vStartTime3<-vStartTime2[vStartTime2[,1]==1,]
  return(list(populationtotal,length(populationtotal),as.vector(vStartTime3[,2])))
}

enrollratef<-c(0.1,0.5,0.9)

####################
repn<-1000
sim11<-vector("list",length(enrollratef))
for (jj in (1:length(enrollratef))){
  seeds<-1:repn
  sim11[[jj]]<-lapply(1:repn,function(x) {
    set.seed(seeds[x])
    pop(vPatsPerMonth=10,nMaxQtyPats=50000,enrollrate1=enrollratef[jj])})
}#if there are 5000 patients in the trial, how many patients will have disease in the population


set.seed(12345)

#################Table 1
Gittins<-GI::Gittins07

vector11<- simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=2,2),K=2,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=30,ptrue=c(0.1,0.3),block=2,criteria=0,rule='FLGI PM',options=0)
####ENS = 7.359 SE = 2.758211 (in the paper 7.7)
all.equal(vector11[[1]], 7.7,tolerance=3*2.758211/sqrt(1000))

vector21<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=2,2),K=2,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=30,ptrue=c(0.1,0.3),block=1,criteria=0,rule='FLGI PM',options=0)
####ENS = 7.61 SE = 2.873353 (in the paper 7.56)
all.equal(vector21[[1]], 7.56,tolerance=3*2.873353/sqrt(1000))

vector31<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=2,2),K=2,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=30,ptrue=c(0.1,0.1),block=2,criteria=0,rule='FLGI PM',options=0)
####ENS = 3.004 SE = 1.653698 (in the paper 3.06)
all.equal(vector31[[1]], 3.06,tolerance=3*1.653698/sqrt(1000))

vector41<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=2,2),K=2,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=30,ptrue=c(0.1,0.1),block=1,criteria=0,rule='FLGI PM',options=0)
####ENS = 3.049 SE = 1.641125  (in the paper 2.99)
all.equal(vector41[[1]], 2.99,tolerance=3*1.641125/sqrt(1000))

vector51<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=2,2),K=2,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=30,ptrue=c(0.7,0.8),block=2,criteria=0,rule='FLGI PM',options=0)
####ENS =  22.82 SE = 2.645486 (in the paper 22.67)
all.equal(vector51[[1]], 22.67,tolerance=3*2.645486/sqrt(1000))

#################Table 3
set.seed(12345)
Gittins<-GI::Gittins099

vector101<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=417,ptrue=c(0.29,0.29,0.29,0.29),block=9,criteria=1,rule='FLGI PM',options=0)

# ENS = 120.7  (in the paper 120.96)
# SE = 9.285916 (in the paper 9.26)
# alpha =  0.048 two sided test 0.004 one sided (in the paper 0.046)
# original code suggested they use two-sided test, which can be seen
# from the matlab example 2. They use meanoc which was two-sided!
# Their matlab code return 0.0036 for one sided test.
# p_star = 0.2493718 (in the paper 0.251)
# SE = 0.2158216 (in the paper 0.21)
all.equal(vector101[[1]], 120.96,tolerance=3*9.285916/sqrt(1000))
all.equal(vector101[[3]], 0.046,tolerance=3*sqrt(0.048*(1-0.048)/1000),scale=1) 
all.equal(vector101[[4]], 0.251,tolerance=3*0.2158216/sqrt(1000))

set.seed(12345)
vector201<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=417,ptrue=c(0.29,0.458, 0.168, 0.24),block=9,criteria=1,rule='FLGI PM',options=0)
# ENS = 179.73(in the paper 179.64)
# SE = 14.69248 (in the paper 13.7)
# alpha = 0.16 (in the paper 0.178)
# p_star = 0.8454527 (in the paper 0.847)
# SE = 0.1229016 (in the paper 0.11)
all.equal(vector201[[1]], 179.64,tolerance=3*14.69248/sqrt(1000),scale=1)
all.equal(vector201[[3]], 0.178,tolerance=3*sqrt(0.16*(1-0.16)/1000),scale=1) 
all.equal(vector201[[4]], 0.847,tolerance=3*0.1229016/sqrt(1000),scale=1)

set.seed(12345)
vector301<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=417,ptrue=c(0.29,0.29,0.29,0.29),block=9,criteria=0,rule='Controlled FLGI',options=0)
# ENS = 120.942 (in the paper 120.86)
# SE = 9.064149 (in the paper 9.22)
# alpha =  0.055  two sided ; 0.032 one sided (in the paper 0.034)
# original code suggested they use one-sided test, which can be seen
# from the matlab example 2. They used meanoc which was one-sided!
# p_star = 0.2506282 (in the paper 0.25)
# SE = 0.02073045 (in the paper 0.02)
all.equal(vector301[[1]], 120.86,tolerance=3*9.064149 /sqrt(1000),scale=1)
all.equal(vector301[[3]], 0.034,tolerance=3*sqrt(0.055*(1-0.055)/1000),scale=1) 
all.equal(vector301[[4]], 0.25,tolerance=3*0.02073045/sqrt(1000),scale=1)

set.seed(44556)
vector401<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                                   Tsize=417,ptrue=c(0.29,0.458, 0.168, 0.24),block=9,criteria=0,rule='Controlled FLGI',options=0)
# ENS = 166.572 (in the paper 166.4)
# SE = 11.86928 (in the paper 11.9)
# alpha = 0.811 (in the paper 0.82)
# p_star = 0.6649261 (in the paper 0.654)
# SE =  0.05933863 (in the paper 0.06)
# Their matlab code 
all.equal(vector401[[1]], 166.4,tolerance=3*11.9/sqrt(1000),scale=1)
all.equal(vector401[[3]], 0.82,tolerance=3*sqrt(0.843*(1-0.843)/1000),scale=1) 
all.equal(vector401[[4]], 0.654,tolerance=3*0.06059153/sqrt(1000),scale=1) 
#p_star not match the table, but I checked their matlab code.
#Their matlab code return:
# ENS = 166.5180
# SE = 11.6284
# alpha = 0.8144
# p_star = 0.6665
# SE =  0.0572

set.seed(12345)
vector501<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=417,ptrue=c(0.458,0.29,0.168, 0.24),block=9,criteria=0,rule='Controlled FLGI',options=0)
# ENS = 130.341 (in the paper 129.84)
# SE = 10.69781 (in the paper 11.2)
# alpha = 0.812 (in the paper 0.812)
# p_star = 0.2494342 (in the paper 0.25)
# SE = 0.02038994 (in the paper 0.02)
all.equal(vector501[[1]], 129.84,tolerance=3*10.69781/sqrt(1000),scale=1)
all.equal(vector501[[3]], 0.812,tolerance=3*sqrt(0.812*(1-0.812)/1000),scale=1) 
all.equal(vector501[[4]], 0.25,tolerance=3*0.02038994/sqrt(1000),scale=1)

set.seed(12345)
vector601<-simula_ocs_adapt_random(jj=1,I0=matrix(1,nrow=4,2),K=4,crit=0.05,noRuns2=100,noRuns=1000,
                        Tsize=417,ptrue=c(0.168,0.458,0.29,0.24),block=9,criteria=0,rule='Controlled FLGI',options=0)
# ENS = 152.58 (in the paper 152.85)
# SE = 12.31435 (in the paper 12.6)
# alpha = 0.988 (in the paper 0.987)
# p_star = 0.6324342 (in the paper 0.634)
# SE =0.09812354 (in the paper 0.1)
all.equal(vector601[[1]], 152.85,tolerance=3*12.6/sqrt(1000),scale=1)
all.equal(vector601[[3]], 0.987,tolerance=3*sqrt(0.987*(1-0.987)/1000),scale=1) 
all.equal(vector601[[4]], 0.634,tolerance=3*0.09812354/sqrt(1000),scale=1)

