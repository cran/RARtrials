
repn<-10

##########block size =2 ##############
set.seed(987650)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.159,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector1=sum(h0decision)/10
all.equal(vector1, 0,tolerance=1e-6,scale=1) 


simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.159,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector2=sum(h1decision)/10
all.equal(vector2, 0.4,tolerance=1e-6,scale=1) 

##########block size =9 ##############
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.0450,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector3=sum(h0decision)/10
all.equal(vector3, 0.1,tolerance=1e-6,scale=1) 


simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.0450,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector4=sum(h1decision)/10
all.equal(vector4, 0.4,tolerance=1e-6,scale=1) 

##########block size =1 ##############

set.seed(987654)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.182,mean=c(0.155,0.155),sd=c(0.64,0.64), side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector5=sum(h0decision)/10
all.equal(vector5, 0,tolerance=1e-6,scale=1) 


set.seed(987655)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=2.182,mean=c(0.155,0.529),sd=c(0.64,0.64), side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector6=sum(h1decision)/10
all.equal(vector6, 0.1,tolerance=1e-6,scale=1) 


##########block size =18 ##############
set.seed(987656)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=1.898,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector7=sum(h0decision)/10
all.equal(vector7, 0.2,tolerance=1e-6,scale=1) 


set.seed(987657)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=rep(2,2),prior_mean=rep(0,2),prior_sd=rep(1,2),
                       stopbound=1.898,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector8=sum(h1decision)/10
all.equal(vector8, 0.4,tolerance=1e-6,scale=1) 


##########block size =1 ##############

set.seed(987658)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=1,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.731,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346) ,side='upper')

})

h0decision<-t(sapply(simnull, "[[", 1))
vector9=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,,drop=F])/10
all.equal(vector9, 0.1,tolerance=1e-6,scale=1) 


set.seed(12345)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=1,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.731,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector10=sum(h1decision[,2]==1)/10
all.equal(vector10, 0.2,tolerance=1e-6,scale=1) 


#######################
set.seed(98765)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=2,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.718,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector101=sum(h1decision[,2]==1)/10
all.equal(vector101, 0.9,tolerance=1e-6,scale=1) 


##########block size =8 ##############
set.seed(9876510)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=8,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.741,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})

h0decision<-t(sapply(simnull, "[[", 1))
vector102=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10
all.equal(vector102, 0.2,tolerance=1e-6,scale=1) 


set.seed(9876511)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=8,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.741,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector103=sum(h1decision[,2]==1)/10
all.equal(vector103, 0.4,tolerance=1e-6,scale=1) 


##########block size =20 ##############
set.seed(9876512)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=20,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.725,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})
h0decision<-t(sapply(simnull, "[[", 1))
vector104=nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,])/10
all.equal(vector104, 0.2,tolerance=1e-6,scale=1) 


set.seed(9876513)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=20,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.725,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})

h1decision<-t(sapply(simalt, "[[", 1))
vector105<-sum(h1decision[,2]==1)/10
all.equal(vector105, 0.7,tolerance=1e-6,scale=1) 

##########block size =40 ##############
set.seed(9876514)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=40,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.662,mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper' )

})
h0decision<-t(sapply(simnull, "[[", 1))
vector106<-nrow(h0decision[h0decision[,2]==1|h0decision[,1]==1,,drop=F])/10
all.equal(vector106, 0.1,tolerance=1e-6,scale=1) 

set.seed(9876515)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=40,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.662,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})
h1decision<-t(sapply(simalt, "[[", 1))
vector107<-sum(h1decision[,2]==1)/10
all.equal(vector107, 0.6,tolerance=1e-6,scale=1) 

##########block size =60 ##############
set.seed(9876515)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_unknown_var(Gittinstype='UNKV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=3,noRuns2=100,
                       Tsize=120,block=60,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),prior_sd=rep(1,3),
                       stopbound=1.591,mean=c(-0.05,0.07,0.13),sd=c(0.346,0.346,0.346) ,side='upper')

})
h1decision<-t(sapply(simalt, "[[", 1))
vector108<-sum(h1decision[,2]==1)/10
all.equal(vector108, 0.6,tolerance=1e-6,scale=1) 



