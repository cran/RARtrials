repn<-10000
##########FLGI-known Table 1 in Williamson, S. F., & Villar, S. S. (2020). A response-adaptive randomization
##########procedure for multi-armed clinical trials with normally distributed outcomes.
##########Biometrics, 76(1), 197–209. https://doi.org/10.1111/biom.13119

##########block size =2 ##############

set.seed(987650)
simnull<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.969,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector1<-sum(h0decision)/10000
#alpha = 0.0492 (in the paper 0.0505)
all.equal(vector1, 0.0505,tolerance=3*sqrt(vector1*(1-vector1)/10000),scale=1) 



set.seed(987651)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=2,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.969,mean=c(0.155,0.529),sd=c(0.64,0.64),side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector2<-sum(h1decision)/10000
#power = 0.2656 (in the paper 0.2701)
all.equal(vector2, 0.2701,tolerance=3*sqrt(vector2*(1-vector2)/10000),scale=1) 

##########block size =9 ##############
set.seed(987652)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.864,mean=c(0.155,0.155),sd=c(0.64,0.64),side='upper' )
})
h0decision<-sapply(simnull, "[[", 1)
vector3<-sum(h0decision)/10000
#alpha = 0.0496 (in the paper 0.0492)
all.equal(vector3, 0.0492,tolerance=3*sqrt(vector3*(1-vector3)/10000),scale=1) 


set.seed(987653)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=9,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.864,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector4<-sum(h1decision)/10000
#power = 0.4205 (in the paper 0.4235)
all.equal(vector4, 0.4235,tolerance=3*sqrt(vector4*(1-vector4)/10000),scale=1) 

##########block size =1 ##############

set.seed(987654)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.991,mean=c(0.155,0.155),sd=c(0.64,0.64) ,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector5<-sum(h0decision)/10000
#alpha = 0.0513 (in the paper 0.0526)
all.equal(vector5, 0.0526,tolerance=3*sqrt(vector5*(1-vector5)/10000),scale=1) 


set.seed(987655)
simalt<-lapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=1,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.991,mean=c(0.155,0.529),sd=c(0.64,0.64) ,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector6<-sum(h1decision)/10000
#power = 0.2289 (in the paper 0.2290)
all.equal(vector6, 0.2290,tolerance=3*sqrt(vector6*(1-vector6)/10000),scale=1) 

##########block size =18 ##############
set.seed(987656)
simnull<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.766,mean=c(0.155,0.155),sd=c(0.64,0.64) ,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector7<-sum(h0decision)/10000
#alpha = 0.0527  (in the paper 0.0513)
all.equal(vector7, 0.0513,tolerance=3*sqrt(vector7*(1-vector7)/10000),scale=1) 


set.seed(987657)
simalt<-parallel::mclapply(1:repn,function(x) {
  sim_flgi_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
                       Tsize=72,block=18,rule='FLGI PM',prior_n=c(1,1),prior_mean=c(0,0),
                       stopbound=1.766,mean=c(0.155,0.529),sd=c(0.64,0.64),side='upper' )
})
h1decision<-sapply(simalt, "[[", 1)
vector8<-sum(h1decision)/10000
#power = 0.5699 (in the paper 0.5653)
all.equal(vector8, 0.5653,tolerance=3*sqrt(vector8*(1-vector8)/10000),scale=1) 

