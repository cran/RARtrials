library(RARfreq)

dabcd_small_var_test<-function(NN,Ntotal1,armn,type,dabcd=0,gamma=2){

  NN<-NN+1/armn
  Ntotal11<-Ntotal1+1#+2
  p<-cbind(p=unname(unlist(NN/Ntotal11)),arm=1:armn)

  rho<-rep(NA,armn)
  rho1<-rep(NA,armn)

  if ( type=='Neyman'){
    for (k in 1:armn ){
      rho1[k]<-sqrt(p[k,1]*(1-p[k,1]))
    }
    for (k in 1:armn ){
      rho[k]<-rho1[k]/(sum(rho1))
    }
  }
  if (type=='RSIHR'){
    for (k in 1:armn ){
      rho1[k]<-sqrt(p[k,1])
    }
    for (k in 1:armn ){
      rho[k]<-rho1[k]/(sum(rho1))
    }
  }

  if (dabcd==1){
    alr<-rep(NA,armn)
    phi<-rep(NA,armn)


    for (k in 1:armn){
      phi[k]<-rho[k]*((rho[k]/(Ntotal1[k]/(sum(Ntotal1)-1)))^gamma)
    }
    for (kk in 1:armn){
      alr[kk]<-phi[kk]/sum(phi)
    }
    return(alr)
  }else if (dabcd==0){
    return(rho)
  }


}




vector11=DBCD_BINARY(S_RK = c(0.4, 0.28, 0.6), N_RK = c(25, 25, 25), rho_func_index = 3, alpha=2)#0.2936500 0.1755583 0.5307917
vector21=dabcd_small_var_test(NN=c(25, 25, 25)*c(0.4, 0.28, 0.6),Ntotal1=c(25, 25, 25),armn=3,type='RSIHR',dabcd=1,gamma=2)#0.2936500 0.1755583 0.5307917
all.equal(vector11$allo_prob, vector21)

vector31=DBCD_BINARY(S_RK = c(0.39, 0.43, 0.4), N_RK = c(82, 79, 88), rho_func_index = 2, alpha=2)#0.3326259 0.3747244 0.2926497
vector41=dabcd_small_var_test(NN=c(82, 79, 88)*c(0.39, 0.43, 0.4),Ntotal1=c(82, 79, 88), armn=3,type='Neyman',dabcd=1,gamma=2)#0.3326259 0.3747244 0.2926497
all.equal(vector31$allo_prob, vector41)

vector51=DBCD_BINARY(S_RK = c(0.4, 0.3,0.38,0.3,0.36),N_RK = c(34,45,65,45,57),rho_func_index = 3, alpha=2)#0.4324360 0.1620624 0.1107438 0.1620624 0.1326954
vector61=dabcd_small_var_test(NN=c(34,45,65,45,57)*c(0.4, 0.3,0.38,0.3,0.36),Ntotal1=c(34,45,65,45,57),armn=5,type='RSIHR',dabcd=1,gamma=2)#0.4324360 0.1620624 0.1107438 0.1620624 0.1326954
all.equal(vector51$allo_prob, vector61)

vector71=DBCD_BINARY(S_RK = c(0.4, 0.3,0.38,0.3,0.36),N_RK = c(34,45,65,45,57),rho_func_index = 2, alpha=2)#0.3940242 0.1843098 0.1051879 0.1843098 0.1321683
vector81=dabcd_small_var_test(NN=c(34,45,65,45,57)*c(0.4, 0.3,0.38,0.3,0.36),Ntotal1=c(34,45,65,45,57),armn=5,type='Neyman',dabcd=1,gamma=2)#0.3940242 0.1843098 0.1051879 0.1843098 0.1321683
all.equal(vector71$allo_prob, vector81)

vector91=DBCD_BINARY(S_RK = c(0.8, 0.7,0.67,0.78,0.75),N_RK = c(77,80,83,90,94),rho_func_index = 2, alpha=2)#0.1867491 0.2537555 0.2535042 0.1500787 0.1559124
vector101=dabcd_small_var_test(NN=c(77,80,83,90,94)* c(0.8, 0.7,0.67,0.78,0.75),Ntotal1=c(77,80,83,90,94),armn=5,type='Neyman',dabcd=1,gamma=2)#0.1867491 0.2537555 0.2535042 0.1500787 0.1559124
all.equal(vector91$allo_prob, vector101)

vector201=DBCD_BINARY(S_RK = c(0.8, 0.7,0.67,0.78,0.75),N_RK = c(77,80,83,90,94),rho_func_index = 3, alpha=2)#0.2675839 0.2031410 0.1768453 0.1889785 0.1634514
vector301=dabcd_small_var_test(NN=c(77,80,83,90,94)* c(0.8, 0.7,0.67,0.78,0.75),Ntotal1=c(77,80,83,90,94),armn=5,type='RSIHR',dabcd=1,gamma=2)#0.2675839 0.2031410 0.1768453 0.1889785 0.1634514
all.equal(vector201$allo_prob, vector301)