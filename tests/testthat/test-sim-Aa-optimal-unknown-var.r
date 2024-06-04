
set.seed(12345)

sim1<-lapply(1:50,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N1=88,N2=880,armn=2,mean=c(9.1/50,8.47/50),sd=c(0.009,0.009),alphaa=0.025,side='lower',
armlabel = c(1,2))})
sim11<-do.call(rbind,lapply(1:50,function(x) {sim1[[x]][[3]]}))
vector11<-table(sim11[,4])/44000
all.equal(as.vector(vector11), c(0.4938636,0.5061364),tolerance=1e-6,scale=1)

set.seed(123451)
sim2<-lapply(1:50,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N1=51,N2=500,armn=3,mean=c(9.1/50,8.47/50,8.8/50),sd=c(0.009,0.01,0.007),alphaa=0.025,side='lower',
armlabel = c(1,2,3))})
sim22<-do.call(rbind,lapply(1:50,function(x) {sim2[[x]][[3]]}))
vector21<-table(sim22[,4])/25000
all.equal(as.vector(vector21), c(0.41680,0.34076,0.24244 ),tolerance=1e-6,scale=1)

set.seed(123452)
sim3<-lapply(1:50,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N1=40,N2=400,armn=4,mean=c(9.1/50,8.47/50,8.8/50,8.6/50),sd=c(0.009,0.01,0.007,0.008),alphaa=0.025,side='upper',
armlabel = c(1,2,3,4))})
sim33<-do.call(rbind,lapply(1:50,function(x) {sim3[[x]][[3]]}))
vector31<-table(sim33[,4])/20000
all.equal(as.vector(vector31), c(0.36745,0.25230,0.17070,0.20955 ),tolerance=1e-6,scale=1)

set.seed(123453)
sim4<-lapply(1:50,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N1=40,N2=400,armn=5,mean=c(9.1/50,8.47/50,8.8/50,8.6/50,8.8/50),sd=c(0.009,0.01,0.007,0.008,0.01),alphaa=0.025,side='upper',
armlabel = c(1,2,3,4,5))})
sim44<-do.call(rbind,lapply(1:50,function(x) {sim4[[x]][[3]]}))
vector41<-table(sim44[,4])/20000
all.equal(as.vector(vector41), c(0.33660,0.18555,0.13465,0.15725,0.18595),tolerance=1e-6,scale=1)




