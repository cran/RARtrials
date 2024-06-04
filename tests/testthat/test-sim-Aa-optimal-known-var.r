
set.seed(12345)

sim1<-lapply(1:100,function(x) {sim_Aa_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N2=880,armn=2,mean=c(9.1/100,8.47/100),sd=c(0.009,0.009),alphaa=0.025,side='lower',armlabel =c(1,2))})
sim11<-do.call(rbind,lapply(1:100,function(x) {sim1[[x]][[3]]}))
vector11<-table(sim11[,4])/88000
all.equal(as.vector(vector11), c(0.4990341,0.5009659),tolerance=1e-6,scale=1)

set.seed(123451)
sim2<-lapply(1:100,function(x) {sim_Aa_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N2=500,armn=3,mean=c(9.1/100,8.47/100,8.8/100),sd=c(0.009,0.01,0.007),alphaa=0.025,side='lower',armlabel =c(1,2,3))})
sim22<-do.call(rbind,lapply(1:100,function(x) {sim2[[x]][[3]]}))
vector21<-table(sim22[,4])/50000
all.equal(as.vector(vector21), c(0.43030,0.33844,0.23126 ),tolerance=1e-6,scale=1)

set.seed(123452)
sim3<-lapply(1:100,function(x) {sim_Aa_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
enrollrate=0.5,N2=400,armn=4,mean=c(9.1/100,8.47/100,8.8/100,8.6/100),sd=c(0.009,0.01,0.007,0.008),alphaa=0.025,side='upper',armlabel =c(1,2,3,4))})
sim33<-do.call(rbind,lapply(1:100,function(x) {sim3[[x]][[3]]}))
vector31<-table(sim33[,4])/40000
all.equal(as.vector(vector31), c(0.384925,0.245675,0.170150,0.199250 ),tolerance=1e-6,scale=1)

set.seed(123453)
sim4<-lapply(1:100,function(x) {sim_Aa_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
enrollrate=0.5,N2=400,armn=5,mean=c(9.1/100,8.47/100,8.8/100,8.6/100,8.8/100),sd=c(0.009,0.01,0.007,0.008,0.01),alphaa=0.025,side='upper',armlabel =c(1,2,3,4,5))})
sim44<-do.call(rbind,lapply(1:100,function(x) {sim4[[x]][[3]]}))
vector41<-table(sim44[,4])/40000
all.equal(as.vector(vector41), c(0.335275,0.188725,0.132225,0.154900,0.188875 ),tolerance=1e-6,scale=1)

