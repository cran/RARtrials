########The set up of the code is try to match some of the results in Table 1, Table 2, Table 5 and Table 6 from
########Guimaraes, P., & Palesch, Y. (2007). Power and sample size simulations for Randomized
########Play-the-Winner rules. Contemporary clinical trials, 28(4), 487–499.
########https://doi.org/10.1016/j.cct.2007.01.006
########Some of the results not match the paper!

PTW<-  function(Pats,nMax,TimeToOutcome,enrollrate,na0,nb0,na1,nb1,ha,hb,Z,N2){#Z2,
  start<-c(rep(1,nb0),rep(0,na0)) #1-treatment;0-control
  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome )
  data1<-matrix(NA_real_,nrow=N2,ncol=5)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime

  phi<-NA

  for (i in 1:N2) {
    data1[i,4]<-sample(start, 1, replace = TRUE)#1-treatment;0-control
    #  data1[i,4]<-sample(start, 1, replace = F)
    if (data1[i,4]==1) {
      data1[i,5]<-rbinom(1,size=1,prob=hb) #1 survival,0 death

      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
      if (total1>0) {
          dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[i,2])),],ncol=5)#data1[1:total1,,drop=F]

        #dataa<-data1[1:total1,,drop=F]
        alive1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==1,,drop=F])#treatment alive
        alive2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==1,,drop=F])#control alive
        deathb<-nrow(dataa[dataa[,4]==1 & dataa[,5]==0,,drop=F])#treatment dead
        death2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==0,,drop=F])#control dead

      }else if (total1==0){
        deathb<-0
        death2<-0
        alive1<-0
        alive2<-0
      }

      start<-c(rep(0,(alive2+deathb)*na1+na0),rep(1,(death2+alive1)*nb1+nb0))

    } else if (data1[i,4]==0) {
      data1[i,5]<-rbinom(1,size=1,prob=ha)
      total1<-sum(as.numeric(data1[,3])<=as.numeric(data1[i,2]))
      if (total1>0) {
        dataa<-matrix(data1[which(as.numeric(data1[,3])<=as.numeric(data1[i,2])),],ncol=5)#data1[1:total1,,drop=F]

        alive1<-nrow(dataa[dataa[,4]==1 & dataa[,5]==1,,drop=F])
        alive2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==1,,drop=F])
        deathb<-nrow(dataa[dataa[,4]==1 & dataa[,5]==0,,drop=F])
        death2<-nrow(dataa[dataa[,4]==0 & dataa[,5]==0,,drop=F])

      }else if (total1==0){
        deathb<-0
        death2<-0
        alive1<-0
        alive2<-0
      }

      start<-c(rep(0,(alive2+deathb)*na1+na0),rep(1,(death2+alive1)*nb1+nb0))
    }
  }
  result<-data.frame(coutcome=data1[,4],doutcome=data1[,5])
  na=nrow(result[ which( result$coutcome==0),])#control
  nb=nrow(result[ which( result$coutcome==1),])
  pa=sum(nrow(result[ which( result$coutcome==0 & result$doutcome==1) , ]))/sum(nrow(result[ which( result$coutcome==0),]))
  pb=sum(nrow(result[ which( result$coutcome==1 & result$doutcome==1) , ]))/sum(nrow(result[ which( result$coutcome==1),]))

  if(is.na(pa)){pa<-0}
  if(is.na(pb)){pb<-0}
  phi<-(pb-pa)^2/(pa*(1-pa)/na+pb*(1-pb)/nb)
   # phi<-(abs(pb-pa)-0.5*(1/na+1/nb))/sqrt(pa*(1-pa)/na+pb*(1-pb)/nb)
  if ( (phi> Z )& !is.nan(phi)) {
    decision<-"Treatment"
  }else {
    decision<-"Control"
  }

  return(list(decision,phi,data1))
}


set.seed(12345)

sim1a<-lapply(1:10000, function(x){
  print(x)
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.1,hb=0.1,Z=2.1488,N2=53)}
)
vector11<-sum(sapply(sim1a, "[[", 2)>(2.143)^2,na.rm=T)/10000
all.equal(vector11, 0.05,tolerance=3*sqrt(vector11*(1-vector11)/10000),scale=1)
#select the Z test statistics to attain type I error around 5% using the sample size around 53
#and in this case the type I error is 0.05


sim1b<-lapply(1:10000, function(x){
  print(x)
   PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
   ha=0.1,hb=0.5,Z=(2.143)^2,N2=53)}
)
vector21<-sum(sapply(sim1b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power  0.9039
#which is close to 0.9
all.equal(vector21, 0.9,tolerance=3*sqrt(vector21*(1-vector21)/10000),scale=1)

set.seed(1234561)
sim2a<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.3,hb=0.3,Z=qnorm(0.975),N2=117)}
)
vector31<-sum(sapply(sim2a, "[[", 2)> (2.015)^2 ,na.rm=T)/10000
#select the Z test statistics to attain type I error around 5% using the sample size around 117
#and in this case the type I error is 0.05
all.equal(vector31, 0.05,tolerance=3*sqrt(vector31*(1-vector31)/10000),scale=1)


sim2b<-lapply(1:10000, function(x){
PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
    ha=0.3,hb=0.6,Z=(2.015)^2,N2=117)})
vector41<-sum(sapply(sim2b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power 0.8948
#which is close to 0.9
all.equal(vector41, 0.9,tolerance=3*sqrt(vector41*(1-vector41)/10000),scale=1)


set.seed(1234562)
sim3a<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.5,hb=0.5,Z=qnorm(0.975),N2=192)}
)
vector51<-sum(sapply(sim3a, "[[", 2)>(1.997)^2 ,na.rm=T)/10000
#select the Z test statistics to attain type I error around 5% using the sample size around 192
#and in this case the type I error is 0.05
all.equal(vector51, 0.05,tolerance=3*sqrt(vector51*(1-vector51)/10000))


sim3b<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.5,hb=0.7,Z=(1.997)^2,N2=192)})
vector61<-sum(sapply(sim3b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power 0.7868
#which is not close to 0.8
all.equal(vector61, 0.8,tolerance=3*sqrt(vector61*(1-vector61)/10000),scale=1)


set.seed(12345)
sim4a<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=2,nb0=2,na1=3,nb1=3,
      ha=0.8,hb=0.8,Z=qnorm(0.95),N2=443)}
)
vector71<-sum(sapply(sim4a, "[[", 2)>(1.947)^2 ,na.rm=T)/10000
#select the Z test statistics to attain type I error around 5% using the sample size around 443
#and in this case the type I error is 0.05
all.equal(vector71, 0.05,tolerance=3*sqrt(vector71*(1-vector71)/10000))


sim4b<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=2,nb0=2,na1=3,nb1=3,
      ha=0.8,hb=0.9,Z=(1.947)^2,N2=443)})
vector81<-sum(sapply(sim4b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power 0.7841
#which is not very close to 0.8
all.equal(vector81, 0.8,tolerance=3*sqrt(vector81*(1-vector81)/10000))



sim5a<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.1,hb=0.1,Z=qnorm(0.975),N2=123)}
)
vector91<-sum(sapply(sim5a, "[[", 2)> (2.03401)^2 ,na.rm=T)/10000
#select the Z test statistics to attain type I error around 5% using the sample size around 443
#and in this case the type I error is 0.0497
all.equal(vector91, 0.05,tolerance=3*sqrt(vector91*(1-vector91)/10000))

sim5b<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.1,hb=0.3,Z=(2.03401)^2,N2=123)})
vector101<-sum(sapply(sim5b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power 0.8084
#which is not close to 0.8
all.equal(vector101, 0.8,tolerance=3*sqrt(vector101*(1-vector101)/10000))


sim6a<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.7,hb=0.7,Z=qnorm(0.975),N2=597)}
)

vector201<-sum(sapply(sim6a, "[[", 2)>(1.9718)^2 ,na.rm=T)/10000
#select the Z test statistics to attain type I error around 5% using the sample size around 597
#and in this case the type I error is 0.05


sim6b<-lapply(1:10000, function(x){
  PTW(Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,na0=1,nb0=1,na1=1,nb1=1,
      ha=0.7,hb=0.8,Z=(1.9718)^2,N2=597)})
vector301<-sum(sapply(sim6b, "[[", 1)=="Treatment")/10000
#using the Z test statistics selected in the previous procedure to attain the power 0.784
#which is not close to 0.8
all.equal(vector301, 0.8,tolerance=3*sqrt(vector301*(1-vector301)/10000),scale=1)
