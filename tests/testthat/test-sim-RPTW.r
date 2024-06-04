set.seed(12345)


sim1<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      h=c(0.3,0.6),Z=qnorm(0.975),N2=117,side='upper')}
)
vector11<-sum(sapply(sim1, "[[", 1)==1)/100
all.equal(vector11, 0.92,tolerance=1e-6,scale=1)


sim2<-lapply(1:100, function(x){
sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
    h=c(0.3,0.6),Z=qnorm(0.975),N2=117,side='upper')})
vector21<-sum(sapply(sim2, "[[", 1)==1)/100
all.equal(vector21, 0.9,tolerance=1e-6,scale=1)


sim3<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      h=c(0.5,0.5),Z=qnorm(0.975),N2=192,side='upper')}
)
vector31<-sum(sapply(sim3, "[[", 1)==1)/100
all.equal(vector31, 0.03,tolerance=1e-6,scale=1)


sim4<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
           h=c(0.5,0.7),Z=qnorm(0.975),N2=192,side='upper')})
vector41<-sum(sapply(sim4, "[[", 1)==1)/100
all.equal(vector41, 0.84,tolerance=1e-6,scale=1)


sim5<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=2,nb0=2,na1=3,nb1=3,
      h=c(0.8,0.8),Z=qnorm(0.975),N2=443,side='upper')}
)
vector51<-sum(sapply(sim5, "[[", 1)==1)/100
all.equal(vector51, 0.04,tolerance=1e-6,scale=1)

sim6<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=2,nb0=2,na1=3,nb1=3,
      h=c(0.8,0.9),Z=qnorm(0.975),N2=443,side='upper')})
vector61<-sum(sapply(sim6, "[[", 1)==1)/100
all.equal(vector61, 0.78,tolerance=1e-6,scale=1)


sim7<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
           h=c(0.1,0.1),Z=qnorm(0.975),N2=123,side='upper')}
)
vector71<-sum(sapply(sim7, "[[", 1)==1)/100
all.equal(vector71, 0.01,tolerance=1e-6,scale=1)

sim8<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
           h=c(0.1,0.3),Z=qnorm(0.975),N2=123,side='upper')})
vector81<-sum(sapply(sim8, "[[", 1)==1)/100
all.equal(vector81, 0.89,tolerance=1e-6,scale=1)

sim9<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
           h=c(0.7,0.7),Z=qnorm(0.975),N2=597,side='upper')}
)
vector91<-sum(sapply(sim9, "[[", 1)==1)/100
all.equal(vector91, 0.02,tolerance=1e-6,scale=1)


sim10<-lapply(1:100, function(x){
  sim_RPTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
           h=c(0.7,0.8),Z=qnorm(0.975),N2=597,side='upper')})
vector101<-sum(sapply(sim10, "[[", 1)==1)/100
all.equal(vector101, 0.78,tolerance=1e-6,scale=1)
