set.seed(12345)

sim1<-lapply(1:1000,function(x) {sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N2=400,armn=2,mean=c(9.1/100,8.47/100),sd=c(0.009,0.009),alphaa=0.025,
cc=mean(c(9.1/100,8.47/100)), armlabel = c(1,2),N1=40,side='lower')})#The optimal allocation ratio is 0.4302534  0.5697466
sim11<-do.call(rbind,lapply(1:1000,function(x) {sim1[[x]][[3]]}))
vector11<-table(sim11[,4])/400000#0.43644 0.56356
all.equal(as.vector(vector11), c(0.4302534,0.5697466),tolerance=1e-1)

set.seed(123451)
sim2<-lapply(1:1000,function(x) {sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N2=300,armn=3,mean=c(9.1/100,8.47/100,8.8/100),sd=c(0.009,0.007,0.01),alphaa=0.025,
cc=mean(c(9.1/100,8.47/100,8.8/100)), armlabel = c(1,2,3),N1=30,side='lower')})#The optimal allocation ratio is 0.3816316 0.2033561 0.4150124
sim22<-do.call(rbind,lapply(1:1000,function(x) {sim2[[x]][[3]]}))
vector21<-table(sim22[,4])/300000
all.equal(as.vector(vector21), c(0.3816316,0.2033561,0.4150124),tolerance=1e-1)

sim21<-lapply(1:1000,function(x) {sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.5,N2=300,armn=3,mean=c(9.1/100,8.47/100,8.8/100),sd=c(0.009,0.007,0.01),alphaa=0.025,
cc=mean(c(9.1/100,8.47/100,8.8/100)), armlabel = c(1,2,3),N1=9,side='lower')})
#The optimal allocation ratio is 0.3816316 0.2033561 0.4150124
sim221<-do.call(rbind,lapply(1:1000,function(x) {sim21[[x]][[3]]}))
vector31<-table(sim221[,4])/300000#0.3806733 0.2074600 0.4118667 (closer results as the initialization period is shorter)
all.equal(as.vector(vector31), c(0.3816316,0.2033561,0.4150124),tolerance=1e-2)

set.seed(123452)
sim3<-lapply(1:1000,function(x) {sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.1,N2=411,armn=4,mean=c(9/100,8.9/100,8.76/100,8.67/100),sd=c(0.009,0.0078,0.01,0.0086),alphaa=0.025,
cc=mean(c(9/100,8.9/100,8.76/100,8.67/100)), armlabel = c(1,2,3,4),N1=40,side='lower')})#The optimal allocation ratio is 0.3477986 0.1689946  0.2777689 0.2054379
sim33<-do.call(rbind,lapply(1:1000,function(x) {sim3[[x]][[3]]}))
vector41<-table(sim33[,4])/411000#0.3388905 0.1752725 0.2757859 0.2100511
all.equal(as.vector(vector41), c(0.3477986,0.1689946,0.2777689,0.2054379),tolerance=1e-1)

sim31<-lapply(1:1000,function(x) {sim_RSIHR_optimal_known_var(Pats=10,nMax=50000,TimeToOutcome=0,
enrollrate=0.1,N2=411,armn=4,mean=c(9/100,8.9/100,8.76/100,8.67/100),
sd=c(0.009,0.0078,0.01,0.0086),alphaa=0.025,
cc=mean(c(9/100,8.9/100,8.76/100,8.67/100)), armlabel = c(1,2,3,4),N1=8,side='upper')})#The optimal allocation ratio is 0.3477986 0.1689946  0.2777689 0.2054379
sim331<-do.call(rbind,lapply(1:1000,function(x) {sim31[[x]][[3]]}))
vector51<-table(sim331[,4])/411000#0.3465085 0.1695718 0.2783601 0.2055596  (closer results as the initialization period is shorter)
all.equal(as.vector(vector51), c(0.3477986,0.1689946,0.2777689,0.2054379),tolerance=1e-2)
