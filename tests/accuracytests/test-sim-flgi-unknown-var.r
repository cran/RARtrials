
repn<-10000
##########Table 1 in Williamson, S. F., & Villar, S. S. (2020). A response-adaptive randomization
##########procedure for multi-armed clinical trials with normally distributed outcomes.
##########Biometrics, 76(1), 197–209. https://doi.org/10.1111/biom.13119

##########block size =2 ##############
set.seed(987650)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.159,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector1=sum(h0decision)/10000
#alpha = 0.0477 (in the paper 0.0497)
all.equal(vector1, 0.0497,tolerance=3*sqrt(vector1*(1-vector1)/10000),scale=1) 


simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.159,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector2=sum(h1decision)/10000
#power = 0.3395 (in the paper 0.3432)
all.equal(vector2, 0.3432,tolerance=3*sqrt(vector2*(1-vector2)/10000),scale=1) 

##########block size =9 ##############
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.0450,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector3=sum(h0decision)/10000
#alpha = 0.0512 (in the paper 0.0514)
all.equal(vector3, 0.0514,tolerance=3*sqrt(vector3*(1-vector3)/10000),scale=1) 


simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.0450,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector4=sum(h1decision)/10000
#power = 0.4253 (in the paper 0.4236)
all.equal(vector4, 0.4236,tolerance=3*sqrt(vector4*(1-vector4)/10000),scale=1) 

##########block size =1 ##############

set.seed(987654)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.182,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector5=sum(h0decision)/10000
#alpha = 0.0491 (in the paper 0.0525)
all.equal(vector5, 0.0525,tolerance=3*sqrt(vector5 *(1-vector5)/10000),scale=1) 


set.seed(987655)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.182,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector6=sum(h1decision)/10000
#power = 0.3275 (in the paper 0.3289)
all.equal(vector6, 0.3289,tolerance=3*sqrt(vector6 *(1-vector6 )/10000),scale=1) 


##########block size =18 ##############


set.seed(987656)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=1.898,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector7=sum(h0decision)/10000
#alpha = 0.0531 (in the paper 0.0517)
all.equal(vector7, 0.0517,tolerance=3*sqrt(vector7 *(1-vector7 )/10000),scale=1) 

RNGkind("L'Ecuyer-CMRG")
set.seed(987657)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=1.898,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector8=sum(h1decision)/10000
#power = 0.5347 (in the paper 0.5277)
all.equal(vector8, 0.5277,tolerance=3*sqrt(vector8 *(1-vector8 )/10000),scale=1) 



##########Table S4 in Williamson, S. F., & Villar, S. S. (2020). A response-adaptive randomization
##########procedure for multi-armed clinical trials with normally distributed outcomes.
##########Biometrics, 76(1), 197–209. https://doi.org/10.1111/biom.13119

######### The power result does not match for this part!!!

##########block size =1 ##############

set.seed(987658)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=1,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.731,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346) ,side='upper')

})

h0decision<-t(sapply(simnull, "[[", 1))
vector9=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10000
#alpha = 0.1056  (in the paper 0.0992)
all.equal(vector9, 0.0992,tolerance=3*sqrt( vector9 *(1- vector9 )/10000),scale=1) 


set.seed(12345)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=1,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.731,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector10=sum(h1decision[,2]==1)/10000
#marginal power =0.533  (in the paper 0.4936)

all.equal(vector10, 0.4936,tolerance=3*sqrt( vector10 *(1- vector10 )/10000),scale=1) 


# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,4]==3)
# }
# mean(aa)    #0.6117917
# 
# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,5])
# }
# mean(aa)    #0.09471732


#######################


#set.seed(12345)
set.seed(98765)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=2,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.718,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector101=sum(h1decision[,2]==1)/10000
#marginal power = 0.5535 (in the paper  0.5121  )
all.equal(vector101, 0.5121,tolerance=3*sqrt( vector101 *(1- vector101 )/10000),scale=1) 

# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,4]==3)
# }
# mean(aa)    # 0.6043792
# 
# 
# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,5])
# }
# mean(aa)  # 85.56%   0.09277538

##########block size =8 ##############

set.seed(9876510)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=8,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.741,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})

h0decision<-t(sapply(simnull, "[[", 1))
vector102=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10000
#alpha = 0.098  (in the paper 0.0954)
all.equal(vector102, 0.0954,tolerance=3*sqrt( vector102 *(1- vector102 )/10000),scale=1) 


set.seed(9876511)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=8,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.741,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector103=sum(h1decision[,2]==1)/10000
#marginal power = 0.5433 (in the paper 0.5235)
all.equal(vector103, 0.5235,tolerance=3*sqrt( vector103 *(1- vector103 )/10000),scale=1) 



##########block size =20 ##############

set.seed(9876512)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=20,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.725,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})
h0decision<-t(sapply(simnull, "[[", 1))
vector104=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10000
#alpha = 0.0991  (in the paper 0.1001)
all.equal(vector104, 0.1001,tolerance=3*sqrt( vector104 *(1- vector104 )/10000),scale=1) 


set.seed(9876513)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=20,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.725,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector105<-sum(h1decision[,2]==1)/10000
#marginal power =0.5823 (in the paper 0.5608)
all.equal(vector105, 0.5608,tolerance=3*sqrt( vector105 *(1- vector105 )/10000),scale=1) 

# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,4]==3)
# }
# mean(aa)    #0.5693358
# 
# 
# 
# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,5])
# }
# mean(aa)    #0.08868439


##########block size =40 ##############
set.seed(9876514)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=40,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.662,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})
h0decision<-t(sapply(simnull, "[[", 1))
vector106<-nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10000
#alpha = 0.1025  (in the paper 0.1009)
all.equal(vector106, 0.1009,tolerance=3*sqrt( vector106 *(1- vector106 )/10000),scale=1) 


set.seed(9876515)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=40,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.662,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector107<-sum(h1decision[,2]==1)/10000
#marginal power = 0.6231 (in the paper 0.61)
all.equal(vector107, 0.61,tolerance=3*sqrt( vector107 *(1- vector107 )/10000),scale=1) 



##########block size =60 ##############

set.seed(9876515)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=60,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.591,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector108<-sum(h1decision[,2]==1)/10000
#marginal power =0.6679  (in the paper 0.6697)
all.equal(vector108, 0.6697,tolerance=3*sqrt( vector108 *(1- vector108)/10000),scale=1) 

# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,4]==3)
# }
# mean(aa)    #0.4827925
# 
# 
# aa<-rep(NA_real_,repn)
# for (i in 1:repn){
#   aa[i]<-mean(simalt[[i]][[4]][,5])
# }
# mean(aa)    #0.07552291


#############Additional code validate the type I error and power under the null with block size =1
#############which match the results above.

GI_Normal_unknown<-Gittins(Gittinstype='UNKV',df=0.995)
valifun<-function(x,mean,sd){
  prior_n1<- prior_n
  posterior_mean<-rep(NA,3)
  posterior_sd<-rep(NA,3)
  xxx<-matrix(NA,nrow=120,ncol=2)
  GI<-rep(NA,3)
  meanhat<-rep(NA,3)
  sdhat<-rep(NA,3)
  nn<-rep(NA,3)
  n<-rep(NA,3)

  for (k in 1:3){
    GI_Std  <- GI_Normal_unknown[prior_n[k]]
    GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std
  }

  optimal_action <-sample(c(1,2,3), size = 1, prob =  c(1/3,1/3,1/3)) #which.is.max(GI)


  xxx[1,1]<-optimal_action
  xxx[1,2]<-rnorm(1,mean=mean[optimal_action],sd=sd[optimal_action])


  for (i in 1:119){

    xxx1<-matrix(xxx[!is.na(xxx)],ncol=2)
    for (k in 1:3){
      n[k]<- nrow(xxx1[xxx1[1:i,1]==k ,,drop=F])
    }

    prior_mean<-rep(0,3)
    prior_sd<-rep(1,3)

    for (k in 1:3){

      if (n[k]==0){
        GI_Std  <- GI_Normal_unknown[prior_n[k]]
        GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std


      }else{

        dataa1<-xxx1[xxx1[,1]==k,2,drop=F]
        kn<-nrow(dataa1)-1
        posterior_mean<-rep(NA,kn+1)
        posterior_sd<-rep(NA,kn+1)


        for (kk in 0:kn){

          posterior_mean[kk+1] <- ((prior_n1[k]+kk)*prior_mean[k]+dataa1[kk+1])/(prior_n1[k]+kk+1)
          posterior_sd[kk+1] <- sqrt(((prior_sd[k])^2)*(prior_n1[k]+kk-1)/(prior_n1[k]+kk) +
                                       (dataa1[kk+1]-prior_mean[k])^2/(prior_n1[k]+kk+1))

          prior_n[k]<-prior_n1[k]+kk
          prior_mean[k]<-posterior_mean[kk+1]
          prior_sd[k]<-posterior_sd[kk+1]

        }

        GI_Std  <- GI_Normal_unknown[ prior_n[k]+1]
        GI[k]   <- prior_mean[k] + prior_sd[k]*GI_Std

      }

    }


    optimal_action <-which.is.max(GI)

    xxx[i+1,1]<- optimal_action
    xxx[i+1,2]<-  rnorm(1,mean[optimal_action],sd[optimal_action])
  }

  for (k in 1:3){
    nn[k]=nrow(xxx[xxx[,1]==k,,drop=F])
    dataaa<-matrix(xxx,ncol=2)
    dataa1<-dataaa[dataaa[,1]==k,2]
    meanhat[k]=mean(dataa1)
    sdhat[k]=sd(dataa1)
  }


  zs1=(meanhat[2]-meanhat[1])/sqrt((sdhat[1])^2/nn[1]+(sdhat[2])^2/nn[2])
  zs2=(meanhat[3]-meanhat[1])/sqrt((sdhat[1])^2/nn[1]+(sdhat[3])^2/nn[3])
  if (is.na(zs1)){zs1=0}
  if (is.na(zs2)){zs2=0}
  if(zs1>=1.731 ){
    ap=1
  }else{
    ap=0
  }
  if(zs2>=1.731 ){
    ap1=1
  }else{
    ap1=0
  }
  if (ap1==1|ap==1){
    app=1
  }else{
    app=0
  }
  return(list(ap,ap1,app))
}

set.seed(12345)
repn=10000
prior_mean<-c(0,0,0)
prior_n<-c(2,2,2)
prior_sd<-c(1,1,1)

simnull<-lapply(1:repn,function(x) {
  print(x)
  valifun(x,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346))
})

h0decision<-t(sapply(simnull, "[[", 3))
vector109<-sum(h0decision)/10000
#0.1028
all.equal(vector109, 0.0992,tolerance=3*sqrt( vector109 *(1- vector109 )/10000),scale=1) 

simalt<-lapply(1:repn,function(x) {
  print(x)
  valifun(x,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346))
})

h1decision<-t(sapply(simalt, "[[", 2))
vector110<-sum(h1decision)/10000
#0.527
all.equal(vector110, 0.4936,tolerance=3*sqrt( vector110 *(1- vector110 )/10000),scale=1) 
