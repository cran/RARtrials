
set.seed(12345)
repn=5

########################
simnull<-lapply(1:repn,function(x) {
sim_flgi_binary(Gittinstype='Binary',df=0.5,Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
                enrollrate=0.5,I0=matrix(1,nrow=2,2),K=2,noRuns2=100,
                Tsize=992,ptrue=c(0.6,0.6),block=20,rule='FLGI PM',ztype='unpooled',stopbound=2.019,,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector11<-sum(h0decision)/5
all.equal(vector11, 0,tolerance=1e-6,scale=1)

simalt<-lapply(1:repn,function(x) {
  sim_flgi_binary(Gittinstype='Binary',df=0.5,Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),
                  enrollrate=0.5,I0=matrix(1,nrow=2,2),K=2,noRuns2=100,
                  Tsize=992,ptrue=c(0.6,0.7),block=20,rule='FLGI PM',ztype='unpooled',stopbound=2.019,,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector21<-sum(h1decision)/5
all.equal(vector21, 0.4,tolerance=1e-6,scale=1)


########################
simnull<-lapply(1:repn,function(x) {
  sim_flgi_binary(Gittinstype='Binary',df=0.995,Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),60, 3)),
                  enrollrate=0.9,I0=matrix(1,nrow=3,2),K=3,noRuns2=100,
                  Tsize=1728,ptrue=c(0.6,0.6,0.6),block=20,rule='FLGI PM',ztype='unpooled',stopbound=2.164,side='upper')
})
h0decision<-sapply(simnull, "[[", 1)
vector31<-sum(h0decision)/5
all.equal(vector31, 0,tolerance=1e-6,scale=1)

simalt<-lapply(1:repn,function(x) {
  sim_flgi_binary(Gittinstype='Binary',df=0.995,Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),60, 3)),
                  enrollrate=0.9,I0=matrix(1,nrow=3,2),K=3,noRuns2=100,
                  Tsize=1728,ptrue=c(0.6,0.7,0.6),block=20,rule='FLGI PM',ztype='unpooled',stopbound=2.164,side='upper')
})
h1decision<-sapply(simalt, "[[", 1)
vector41<-sum(h1decision)/5
all.equal(vector41, 0.8,tolerance=1e-6,scale=1)


simalt<-lapply(1:repn,function(x) {
  sim_flgi_binary(Gittinstype='Binary',df=0.995,Pats=10,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),60, 3)),
                  enrollrate=0.9,I0=matrix(1,nrow=3,2),K=3,noRuns2=100,
                  Tsize=1728,ptrue=c(0.6,0.7,0.7),block=20,rule='FLGI PM',ztype='unpooled',stopbound=2.164,side='upper')
})
h1decision<-t(sapply(simalt, "[[", 1))
vector51<-sum(h1decision[,1])/5
vector61<-sum(h1decision[,2])/5
all.equal(vector51, 0.6,tolerance=1e-6,scale=1)
all.equal(vector61, 0.8,tolerance=1e-6,scale=1)

