repn<-10

##########block size =2 ##############
set.seed(987650)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.969,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector1<-sum(h0decision)/10
all.equal(vector1, 0.1,tolerance=1e-6,scale=1)  



set.seed(987651)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.969,mean=c(0.155,0.529),sd=c(0.64,0.64),side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector2<-sum(h1decision)/10
all.equal(vector2, 0.6,tolerance=1e-6,scale=1)  

##########block size =9 ##############
set.seed(987652)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.864,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector3<-sum(h0decision)/10
all.equal(vector3, 0,tolerance=1e-6,scale=1) 


set.seed(987653)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.864,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector4<-sum(h1decision)/10
all.equal(vector4, 0.3,tolerance=1e-6,scale=1)  

##########block size =1 ##############

set.seed(987654)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.991,mean=c(0.155,0.155),sd=c(0.64,0.64) ,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector5<-sum(h0decision)/10
all.equal(vector5, 0,tolerance=1e-6,scale=1) 


set.seed(987655)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.991,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector6<-sum(h1decision)/10
all.equal(vector6, 0.4,tolerance=1e-6,scale=1) 

##########block size =18 ##############
set.seed(987656)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.766,mean=c(0.155,0.155),sd=c(0.64,0.64) ,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector7<-sum(h0decision)/10
all.equal(vector7, 0.1,tolerance=1e-6,scale=1) 


set.seed(987657)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.766,mean=c(0.155,0.529),sd=c(0.64,0.64),side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector8<-sum(h1decision)/10
all.equal(vector8, 0.3,tolerance=1e-6,scale=1) 

