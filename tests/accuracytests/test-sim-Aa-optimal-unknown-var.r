
set.seed(12345)

sim1<-lapply(1:1000,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N1=88,N2=880,armn=2,mean=c(9.1/100,8.47/100),sd=c(0.009,0.009),alphaa=0.025,side='lower',
armlabel = c(1,2))})
sim11<-do.call(rbind,lapply(1:1000,function(x) {sim1[[x]][[3]]}))
vector11<-table(sim11[,4])/880000#allocation ratio should be around 0.5#0.4999227 0.5000773
all.equal(as.vector(vector11), c(0.5,0.5),tolerance=1e-3)

set.seed(123451)
sim2<-lapply(1:1000,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.9,N1=51,N2=500,armn=3,mean=c(9.1/100,8.47/100,8.8/100),sd=c(0.009,0.01,0.007),alphaa=0.025,side='lower',
armlabel = c(1,2,3))})
sim22<-do.call(rbind,lapply(1:1000,function(x) {sim2[[x]][[3]]}))
vector21<-table(sim22[,4])/500000#allocation ratio should be around  0.428147, 0.3363841, 0.2354689 #0.418940 0.335374 0.245686
all.equal(as.vector(vector21), c(0.428147, 0.3363841, 0.2354689),tolerance=1e-1) #2e-2=0.02

set.seed(123452)
sim3<-lapply(1:1000,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N1=40,N2=400,armn=4,mean=c(9.1/100,8.47/100,8.8/100,8.6/100),sd=c(0.009,0.01,0.007,0.008),alphaa=0.025,side='upper',
armlabel = c(1,2,3,4))})
sim33<-do.call(rbind,lapply(1:1000,function(x) {sim3[[x]][[3]]}))
vector31<-table(sim33[,4])/400000#allocation ratio should be around  0.3840613, 0.2463755, 0.1724628, 0.1971004#0.3721475 0.2465225 0.1786800 0.2026500
all.equal(as.vector(vector31), c(0.3840613, 0.2463755, 0.1724628, 0.1971004),tolerance=1e-1)

set.seed(123453)
sim4<-lapply(1:1000,function(x) {sim_Aa_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N1=40,N2=400,armn=5,mean=c(9.1/100,8.47/100,8.8/100,8.6/100,8.8/100),sd=c(0.009,0.01,0.007,0.008,0.01),alphaa=0.025,side='upper',
armlabel = c(1,2,3,4,5))})
sim44<-do.call(rbind,lapply(1:1000,function(x) {sim4[[x]][[3]]}))
vector41<-table(sim44[,4])/400000#allocation ratio should be around  0.3396226415 0.1886792 0.13207547 0.150943396 0.188679#0.3267500 0.1894800 0.1382625 0.1551450 0.1903625
all.equal(as.vector(vector41), c(0.3396226415, 0.1886792, 0.13207547, 0.150943396, 0.188679),tolerance=1e-1)




