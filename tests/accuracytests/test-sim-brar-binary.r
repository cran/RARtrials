set.seed(12345)
repn<-10000

sim<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),1, 0.3)),
                   enrollrate=1,N1=50,armn=5,h=c(0.2,0.2,0.2,0.2,0.4),au=c(0.345,0.345,0.345,0.345),
                   N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                   alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                   alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),tpp=0,
                   deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})


h1decision<-t(sapply(sim, "[[", 1) )

vector11<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",])/10000 #0.0342
vector21<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",])/10000 #0.036
vector31<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",])/10000 #0.0313
vector41<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",])/10000 #0.6686
all.equal(vector11, 0.03,tolerance=sqrt(3)*sqrt(0.03*0.97/10000))
all.equal(vector21, 0.03,tolerance=sqrt(3)*sqrt(0.03*0.97/10000))
all.equal(vector31, 0.03,tolerance=sqrt(3)*sqrt(0.03*0.97/10000))
all.equal(vector41, 0.67,tolerance=sqrt(3)*sqrt(0.03*0.97/10000))

vector51<-nrow(h1decision[h1decision[,2]=="Futility",])/10000 #0.5617
vector61<-nrow(h1decision[h1decision[,3]=="Futility",])/10000 #0.5613
vector71<-nrow(h1decision[h1decision[,4]=="Futility",])/10000 #0.5639
vector81<-nrow(h1decision[h1decision[,5]=="Futility",])/10000 #0.0777
all.equal(vector51, 0.58,tolerance=1e-2)
all.equal(vector61, 0.58,tolerance=1e-2)
all.equal(vector71, 0.58,tolerance=1e-2)
all.equal(vector81, 0.07,tolerance=1e-2)

#####################
#####################
#####################
sim<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),1, 0.3)),
                  enrollrate=1,N1=50,armn=5,h=c(0.2,0.25,0.3,0.35,0.4),au=c(0.345,0.345,0.345,0.345),
                  N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                  alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                  alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),
                  tpp=0,deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})

h1decision<-t(sapply(sim, "[[", 1) )

vector101<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",])/10000 #0.0988
vector201<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",])/10000 #0.2148
vector301<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",])/10000 #0.3732
vector401<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",])/10000 #0.6068
all.equal(vector101, 0.10,tolerance=1e-2)
all.equal(vector201, 0.21,tolerance=1e-2)
all.equal(vector301, 0.38,tolerance=1e-2)
all.equal(vector401, 0.61,tolerance=1e-2)

vector501<-nrow(h1decision[h1decision[,2]=="Futility",])/10000 #0.365
vector601<-nrow(h1decision[h1decision[,3]=="Futility",])/10000 #0.2404
vector701<-nrow(h1decision[h1decision[,4]=="Futility",])/10000 #0.1416
vector801<-nrow(h1decision[h1decision[,5]=="Futility",])/10000 #0.073
all.equal(vector501, 0.36,tolerance=1e-2)
all.equal(vector601, 0.22,tolerance=1e-2)
all.equal(vector701, 0.13,tolerance=1e-2)
all.equal(vector801, 0.06,tolerance=1e-2)

#####################
#####################
#####################
sim<-lapply(1:repn,function(x) {
  sim_brar_binary(Pats=5,nMax=50000,TimeToOutcome=expression(rnorm( length( vStartTime ),1, 0.3)),
                  enrollrate=1,N1=50,armn=5,h=c(0.2,0.2,0.2,0.2,0.2),au=c(0.345,0.345,0.345,0.345),
                  N2=250,tp=1,armlabel=c(1,2,3,4,5),blocksize=10,
                  alpha1=0.2,beta1=0.8,alpha2=0.2,beta2=0.8,alpha3=0.2,beta3=0.8,
                  alpha4=0.2,beta4=0.8,alpha5=0.2,beta5=0.8,minstart=50,deltaa=c(0.2,0.2,0.2,0.2),
                  tpp=0,deltaa1=c(0.2,0.2,0.2,0.2),side='upper')
})


h1decision<-t(sapply(sim, "[[", 1) )

vector1001<-nrow(h1decision[h1decision[,2]=="Superiorityfinal",])/10000 #0.0178
vector2001<-nrow(h1decision[h1decision[,3]=="Superiorityfinal",])/10000 #0.0138
vector3001<-nrow(h1decision[h1decision[,4]=="Superiorityfinal",])/10000 #0.0122
vector4001<-nrow(h1decision[h1decision[,5]=="Superiorityfinal",])/10000 #0.012
all.equal(vector1001, 0.01,tolerance=1e-2)
all.equal(vector2001, 0.01,tolerance=1e-2)
all.equal(vector3001, 0.01,tolerance=1e-2)
all.equal(vector4001, 0.01,tolerance=1e-2)

vector5001<-nrow(h1decision[h1decision[,2]=="Futility",])/10000 #0.7518
vector6001<-nrow(h1decision[h1decision[,3]=="Futility",])/10000 #0.7586
vector7001<-nrow(h1decision[h1decision[,4]=="Futility",])/10000 # 0.7528
vector8001<-nrow(h1decision[h1decision[,5]=="Futility",])/10000 #0.756
all.equal(vector5001, 0.74,tolerance=1e-2)
all.equal(vector6001, 0.74,tolerance=1e-2)
all.equal(vector7001, 0.74,tolerance=1e-2)
all.equal(vector8001, 0.74,tolerance=1e-2)
