set.seed(12345)
repn<-50

sim1<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm(length( vStartTime ),1, 0.3)),
                   enrollrate=1,N1=50,armn=5,h=c(0.2,0.2,0.2,0.2,0.4),au=c(0.345,0.345,0.345,0.345),
                   N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                   alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                   alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),tpp=0,
                   deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})


h1decision<-t(sapply(sim1, "[[", 1) )

vector11<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",,drop=F])/repn 
vector21<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",,drop=F])/repn 
vector31<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",,drop=F])/repn 
vector41<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",,drop=F])/repn 
all.equal(vector11, 0.02,tolerance=1e-6,scale=1)
all.equal(vector21, 0.06,tolerance=1e-6,scale=1)
all.equal(vector31, 0.04,tolerance=1e-6,scale=1)
all.equal(vector41, 0.66,tolerance=1e-6,scale=1)

vector51<-nrow(h1decision[h1decision[,2]=="Futility",,drop=F])/repn 
vector61<-nrow(h1decision[h1decision[,3]=="Futility",,drop=F])/repn 
vector71<-nrow(h1decision[h1decision[,4]=="Futility",,drop=F])/repn 
vector81<-nrow(h1decision[h1decision[,5]=="Futility",,drop=F])/repn 
all.equal(vector51, 0.68,tolerance=1e-6,scale=1)
all.equal(vector61, 0.54,tolerance=1e-6,scale=1)
all.equal(vector71, 0.58,tolerance=1e-6,scale=1)
all.equal(vector81, 0.02,tolerance=1e-6,scale=1)

#####################
#####################
#####################
sim2<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),1, 0.3)),
                  enrollrate=1,N1=50,armn=5,h=c(0.2,0.25,0.3,0.35,0.4),au=c(0.345,0.345,0.345,0.345),
                  N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                  alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                  alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),
                  tpp=0,deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})

h1decision<-t(sapply(sim2, "[[", 1) )

vector101<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",,drop=F])/repn 
vector201<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",,drop=F])/repn 
vector301<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",,drop=F])/repn 
vector401<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",,drop=F])/repn 
all.equal(vector101, 0.08,tolerance=1e-6,scale=1)
all.equal(vector201, 0.34,tolerance=1e-6,scale=1)
all.equal(vector301, 0.36,tolerance=1e-6,scale=1)
all.equal(vector401, 0.64,tolerance=1e-6,scale=1)

vector501<-nrow(h1decision[h1decision[,2]=="Futility",,drop=F])/repn 
vector601<-nrow(h1decision[h1decision[,3]=="Futility",,drop=F])/repn 
vector701<-nrow(h1decision[h1decision[,4]=="Futility",,drop=F])/repn 
vector801<-nrow(h1decision[h1decision[,5]=="Futility",,drop=F])/repn 
all.equal(vector501, 0.26,tolerance=1e-6,scale=1)
all.equal(vector601, 0.28,tolerance=1e-6,scale=1)
all.equal(vector701, 0.24,tolerance=1e-6,scale=1)
all.equal(vector801, 0.06,tolerance=1e-6,scale=1)

#####################
#####################
#####################
sim3<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),1, 0.3)),
                  enrollrate=1,N1=50,armn=5,h=c(0.2,0.2,0.2,0.2,0.2),au=c(0.345,0.345,0.345,0.345),
                  N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                  alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                  alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),
                  tpp=0,deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})


h1decision<-t(sapply(sim3, "[[", 1) )

vector1001<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",,drop=F])/repn 
vector2001<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",,drop=F])/repn 
vector3001<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",,drop=F])/repn 
vector4001<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",,drop=F])/repn 
all.equal(vector1001, 0.02,tolerance=1e-6,scale=1)
all.equal(vector2001, 0.02,tolerance=1e-6,scale=1)
all.equal(vector3001, 0.04,tolerance=1e-6,scale=1)
all.equal(vector4001, 0.02,tolerance=1e-6,scale=1)

vector5001<-nrow(h1decision[h1decision[,2]=="Futility",,drop=F])/repn 
vector6001<-nrow(h1decision[h1decision[,3]=="Futility",,drop=F])/repn 
vector7001<-nrow(h1decision[h1decision[,4]=="Futility",,drop=F])/repn 
vector8001<-nrow(h1decision[h1decision[,5]=="Futility",,drop=F])/repn 
all.equal(vector5001, 0.74,tolerance=1e-6,scale=1)
all.equal(vector6001, 0.78,tolerance=1e-6,scale=1)
all.equal(vector7001, 0.82,tolerance=1e-6,scale=1)
all.equal(vector8001, 0.76,tolerance=1e-6,scale=1)
