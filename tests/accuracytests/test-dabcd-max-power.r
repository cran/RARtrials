####dabcd_max_power_test try to solve the optimizer, but it does not work for discontinuous function and
####only works for several scenarios and is only used for test.
library(nloptr)
dabcd_max_power_test<-function (B,w,p,nn,ub,armn){#

  if (armn==3){
    C=qchisq(0.95,2)
    eval_f <- function(n) {
      return(w[1]*n[1]+w[2]*n[2]+w[3]*n[3])
    }
    eval_g_ineq <- function (n) {

      constr <- c(B*(n[1]+n[2]+n[3])-n[1],
                  B*(n[1]+n[2]+n[3])-n[2],
                  B*(n[1]+n[2]+n[3])-n[3],
                  C+(((n[1]/(p[1]*(1-p[1])))*(p[1]-p[3])+(n[2]/(p[2]*(1-p[2])))*(p[2]-p[3]))^2)/((n[1]/(p[1]*(1-p[1])))+(n[2]/(p[2]*(1-p[2])))+(n[3]/(p[3]*(1-p[3]))))-
                    ((n[1]/(p[1]*(1-p[1])))*(p[1]-p[3])^2)-((n[2]/(p[2]*(1-p[2])))*(p[2]-p[3])^2))
      return (constr)
    }

    opts <- list( "algorithm" = "NLOPT_LN_COBYLA", 
                  "xtol_rel" = 1.0e-15,
                  "maxeval"= 1000,
                  "tol_constraints_ineq" = rep( 1.0e-10, 4 )
                  )

  }else if  (armn==4){
    C=qchisq(0.95,3)
    eval_f <- function(n) {
      return(w[1]*n[1]+w[2]*n[2]+w[3]*n[3]+w[4]*n[4])
    }
    eval_g_ineq <- function (n) {
      constr <- c(B*(n[1]+n[2]+n[3]+n[4])-n[1],
                  B*(n[1]+n[2]+n[3]+n[4])-n[2],
                  B*(n[1]+n[2]+n[3]+n[4])-n[3],
                  B*(n[1]+n[2]+n[3]+n[4])-n[4],
                  C+(((n[1]/(p[1]*(1-p[1])))*(p[1]-p[4])+(n[2]/(p[2]*(1-p[2])))*(p[2]-p[4])+(n[3]/(p[3]*(1-p[3])))*(p[3]-p[4]))^2)/((n[1]/(p[1]*(1-p[1])))+(n[2]/(p[2]*(1-p[2])))+(n[3]/(p[3]*(1-p[3])))+(n[4]/(p[4]*(1-p[4]))))-
                    ((n[1]/(p[1]*(1-p[1])))*(p[1]-p[4])^2)-((n[2]/(p[2]*(1-p[2])))*(p[2]-p[4])^2)-((n[3]/(p[3]*(1-p[3])))*(p[3]-p[4])^2))
      return (constr)
    }

    opts <- list( "algorithm" = "NLOPT_LN_COBYLA", 
                  "xtol_rel" = 1.0e-15,
                  "maxeval"= 1000,
                  "tol_constraints_ineq" = rep( 1.0e-10, 5 ))

  }else if  (armn==5){
    C=qchisq(0.95,4)
    eval_f <- function(n) {
      return(w[1]*n[1]+w[2]*n[2]+w[3]*n[3]+w[4]*n[4]+w[5]*n[5])
    }
    eval_g_ineq <- function (n) {
      constr <- c(B*(n[1]+n[2]+n[3]+n[4]+n[5])-n[1],
                  B*(n[1]+n[2]+n[3]+n[4]+n[5])-n[2],
                  B*(n[1]+n[2]+n[3]+n[4]+n[5])-n[3],
                  B*(n[1]+n[2]+n[3]+n[4]+n[5])-n[4],
                  B*(n[1]+n[2]+n[3]+n[4]+n[5])-n[5],
                  C+(((n[1]/(p[1]*(1-p[1])))*(p[1]-p[5])+(n[2]/(p[2]*(1-p[2])))*(p[2]-p[5])+(n[3]/(p[3]*(1-p[3])))*(p[3]-p[5])+(n[4]/(p[4]*(1-p[4])))*(p[4]-p[5]))^2)/((n[1]/(p[1]*(1-p[1])))+(n[2]/(p[2]*(1-p[2])))+(n[3]/(p[3]*(1-p[3])))+(n[4]/(p[4]*(1-p[4])))+(n[5]/(p[5]*(1-p[5]))))-
                    ((n[1]/(p[1]*(1-p[1])))*(p[1]-p[5])^2)-((n[2]/(p[2]*(1-p[2])))*(p[2]-p[5])^2)-((n[3]/(p[3]*(1-p[3])))*(p[3]-p[5])^2)-((n[4]/(p[4]*(1-p[4])))*(p[4]-p[5])^2))
      return (constr)
    }

    opts <- list( "algorithm" = "NLOPT_LN_COBYLA", 
                  "xtol_rel" = 1.0e-15,
                  "maxeval"= 1000,
                  "tol_constraints_ineq" = rep( 1.0e-10, 6 ))

  }



  x0<-nn #change here
  lb <-rep(1,armn)#c(0,0,0)
  ub <- ub

  res <- nloptr(
    x0 = x0,
    eval_f = eval_f,
   lb = lb,
    ub = ub,
    eval_g_ineq = eval_g_ineq,
    opts = opts )
  pp<-res$solution
  pplist<-vector("list",armn)
  for (j in 1:length(pp)) {
    pplist[[j]]<-pp[[j]]/sum(pp)
  }

  return(as.vector(do.call(rbind,pplist)))
}


#####Three arms Neyman
vector11=dabcd_max_power_test(B=0.25,w=c(1,1,1),p=(c(17,9,19)+1)/(c(129,155,123)+2) ,nn=c(129,155,123)+2,armn=3,ub=c(960,960,960)) #0.2500000 0.3688891 0.3811109
vector21=dabcd_max_power(NN=c(17,9,19),Ntotal1=c(129,155,123),armn=3, BB=0.25, type='Neyman',dabcd=0)#0.2500000 0.3688891 0.3811109
all.equal(vector11, vector21)

vector31=dabcd_max_power_test(B=0.3,w=c(1,1,1),p=(c(1,1,2)+1)/(c(3,5,4)+2) ,nn=c(3,5,4)+2,armn=3,ub=c(960,960,960)) #0.3000000 0.3488257 0.3511743
vector41=dabcd_max_power(NN=c(1,1,2),Ntotal1=c(3,5,4),armn=3, BB=0.3, type='Neyman',dabcd=0) #0.3000000 0.3488257 0.3511743
all.equal(vector31, vector41)

vector51=dabcd_max_power_test(B=0.2,w=c(1,1,1),p=(c(3, 6, 8)+1)/(c(34, 36, 34)+2),nn=c(34, 36, 34)+2,ub=c(230,230,230),armn=3)#0.3555952 0.2000000 0.4444048
vector61=dabcd_max_power(NN=c(3, 6, 8),Ntotal1=c(34, 36, 34),armn=3, BB=0.2,type='Neyman',dabcd=0)#0.3555952 0.2000000 0.4444048
all.equal(vector51, vector61)

vector71=dabcd_max_power_test(B=0.1,w=c(1,1,1),p=(c(30, 20, 18)+1)/(c(90, 100, 105)+2),nn=c(90, 100, 105)+2,ub=c(160,160,160),armn=3)#0.5274177 0.1000000 0.3725823
vector81=dabcd_max_power(NN=c(30, 20, 18),Ntotal1=c(90, 100, 105),armn=3, BB=0.1, type='Neyman',dabcd=0)#0.5274177 0.1000000 0.3725823
all.equal(vector71, vector81)

vector91=dabcd_max_power(NN=c(35,66,50),Ntotal1=c(100,88,90),armn=3,BB=0.1, type='Neyman',dabcd=0)#0.4698286 0.4301714 0.1000000
vector101=dabcd_max_power_test(B=0.1,w=c(1,1,1),p=(c(35,66,50)+1)/(c(100,88,90)+2) ,nn=c(100,88,90)+2,ub=c(460,460,460),armn=3)#0.4698286 0.4301714 0.1000000
all.equal(vector91, vector101)

vector201=dabcd_max_power(NN=c(35,66,66),Ntotal1=c(100,90,90),armn=3,BB=0.1, type='Neyman',dabcd=0) #0.5178969 0.2410515 0.2410515
p=c(0.3529412 ,0.7282609 ,0.7282609)
q=1-p
sqrt(p[1]*q[1])/(sqrt(p[1]*q[1])+sqrt(p[2]*q[2]))#0.517897
sqrt(p[2]*q[2])/2/(sqrt(p[1]*q[1])+sqrt(p[2]*q[2]))#0.2410515
all.equal(vector201, c(0.517897,0.2410515,0.2410515),tolerance=1e-6)

vector301=dabcd_max_power(NN=c(35,66,66),Ntotal1=c(100,90,90),armn=3,BB=0.25, type='Neyman',dabcd=0)#0.50 0.25 0.25
all.equal(vector301, c(0.50,0.25,0.25))

vector401=dabcd_max_power(NN=c(120,120,30),Ntotal1=c(150,150,90),armn=3,BB=0.1, type='Neyman',dabcd=0)# 0.2300881 0.2300881 0.5398239
p=c(0.7960526, 0.7960526, 0.3369565)
q=1-p
sqrt(p[3]*q[3])/(sqrt(p[1]*q[1])+sqrt(p[3]*q[3]))# 0.5398238
sqrt(p[1]*q[1])/2/(sqrt(p[1]*q[1])+sqrt(p[3]*q[3]))# 0.2300881
all.equal(vector401, c(0.2300881,0.2300881,0.5398239),tolerance=1e-6)

vector501=dabcd_max_power(NN=c(120,120,30)+1,Ntotal1=c(150,150,90)+2,armn=3,BB=0.25, type='Neyman',dabcd=0)#0.25 0.25 0.50
all.equal(vector501, c(0.25,0.25,0.50))

#####Three arms RSIHR
vector11=dabcd_max_power(NN=c(40,52,66),Ntotal1=c(102,130,130),armn=3,BB=0.1, type='RSIHR',dabcd=0 )#0.375952 0.100000 0.524048
vector21=dabcd_max_power_test(B=0.1,w=1-(c(40,52,66)+1)/(c(102,130,130)+2) ,p=(c(40,52,66)+1)/(c(102,130,130)+2),
                nn=c(102,130,130)+2,ub=c(250,250,250),armn=3)#0.375952 0.100000 0.524048
all.equal(vector11, vector21)

vector31=dabcd_max_power(NN=c(80,38,55),Ntotal1=c(110,88,96),armn=3,BB=0.2, type='RSIHR',dabcd=0  )# 0.4583773 0.3416227 0.2000000
vector41=dabcd_max_power_test(B=0.2,w=1-(c(80,38,55)+1)/(c(110,88,96)+2) ,p=(c(80,38,55)+1)/(c(110,88,96)+2),
                nn=c(110,88,96)+2,ub=c(180,180,180),armn=3)# 0.4583773 0.3416227 0.2000000
all.equal(vector31, vector41)

vector51=dabcd_max_power(NN=c(80,55,55),Ntotal1=c(110,88,88),armn=3,BB=0.15, type='RSIHR' ,dabcd=0 )#0.5187922 0.2406039 0.2406039
p=c(0.7232143,0.6222222,0.6222222)
sqrt(p[3])/2/(sqrt(p[1])+sqrt(p[3]))#0.2406039
sqrt(p[1])/(sqrt(p[1])+sqrt(p[3]))# 0.5187922
all.equal(vector51, c(0.5187922,0.2406039,0.2406039),tolerance=1e-6)

vector61=dabcd_max_power(NN=c(80,55,55),Ntotal1=c(110,88,88),armn=3,BB=0.25, type='RSIHR')# 0.50 0.25 0.25
all.equal(vector61, c(0.50,0.25,0.25))

#########Four Arms
vector11=dabcd_max_power(NN=c(54,67,85,63),Ntotal1=c(100,88,90,94),armn=4,BB=0.2, type='Neyman')#0.3466227 0.2000000 0.2533773 0.2000000
vector21=dabcd_max_power_test(B=0.2,w=c(1,1,1,1),p=(c(54,67,85,63)+1)/(c(100,88,90,94)+2),nn=c(100,88,90,94)+2,ub=c(500,500,500,500),armn=4)#0.3466227 0.2000000 0.2533773 0.2000000
all.equal(vector11, vector21)

vector31=dabcd_max_power(NN=c(120,80,111,93),Ntotal1=c(154,128,140,160),armn=4,BB=0.1, type='Neyman')#0.1000000 0.1000000 0.3463452 0.4536548
vector41=dabcd_max_power_test(B=0.1,w=c(1,1,1,1),p=(c(120,80,111,93)+1)/(c(154,128,140,160)+2),nn=c(154,128,140,160)+2,ub=c(420,420,420,420),armn=4)#0.1000000 0.1000000 0.3463452 0.4536548
all.equal(vector31, vector41)

vector51=dabcd_max_power(NN=c(80,80,80,93),Ntotal1=c(154,154,154,160),armn=4,BB=0.1, type='Neyman')#0.1676923 0.1676923 0.1676923 0.4969231
p=c(0.5192308,0.5192308,0.5192308,0.5802469)
q=1-p
sqrt(p[4]*q[4])/(sqrt(p[1]*q[1])+sqrt(p[4]*q[4]))# 0.4969231
sqrt(p[1]*q[1])/3/(sqrt(p[1]*q[1])+sqrt(p[4]*q[4]))#  0.1676923
all.equal(vector51, c(0.1676923,0.1676923,0.1676923,0.4969231),tolerance=1e-6)

vector61=dabcd_max_power(NN=c(80,80,80,93),Ntotal1=c(154,154,154,160),armn=4,BB=0.2, type='Neyman')# 0.2 0.2 0.2 0.4
all.equal(vector61, c(0.2,0.2,0.2,0.4),tolerance=1e-6)

vector71=dabcd_max_power(NN=c(54,54,85,63)+1,Ntotal1=c(88,88,90,94)+2,armn=4,BB=0.2, type='Neyman')#0.2326693 0.2326693 0.3346614 0.2000000
vector81=dabcd_max_power_test(B=0.2,w=c(1,1,1),p=as.vector((c(54,85,63)+1)/(c(88,90,94)+2) ),nn=c(88,90,94)+2,ub=c(500,500,500)+2,armn=3)#0.4777227 0.3222773 0.2000000 (0.4777227/2=0.2388614)
all.equal(vector71, vector81)

vector91=dabcd_max_power(NN=c(154,154,185,163),Ntotal1=c(188,188,190,194),armn=4,BB=0.2,type='Neyman')#0.2487669 0.2487669 0.3024662 0.2000000
vector101=dabcd_max_power_test(B=0.2,w=c(1,1,1),p=as.vector((c(154,185,163)+1)/(c(188,190,194)+2) ),nn=c(188,190,194)+2,ub=c(500,500,500)+2,armn=3)#0.4975338 0.3024662 0.2000000 (0.4975338/2=0.2487669)
all.equal(vector91, c(vector101[1]/2,vector101[1]/2,vector101[2],vector101[3]))

vector201=dabcd_max_power(NN=c(120,120,120,30),Ntotal1=c(150,150,150,90),armn=4,BB=0.1, type='Neyman' )#0.1533920 0.1533920 0.1533920 0.5398239
p=c(0.7960526,0.7960526,0.7960526,0.3369565)
q=1-p
sqrt(p[4]*q[4])/(sqrt(p[1]*q[1])+sqrt(p[4]*q[4]))# 0.5398238
sqrt(p[1]*q[1])/3/(sqrt(p[1]*q[1])+sqrt(p[4]*q[4]))#  0.1533921
all.equal(vector201, c(0.1533920,0.1533920,0.1533920,0.5398239),tolerance=1e-6)

vector301=dabcd_max_power(NN=c(120,120,120,30),Ntotal1=c(150,150,150,90),armn=4,BB=0.25, type='Neyman' )#0.25 0.25 0.25 0.25
all.equal(vector301, c(0.25,0.25,0.25,0.25))

vector401=dabcd_max_power(NN=c(85,74,99,99)+1,Ntotal1=c(143,138,160,160)+2,armn=4,BB=0.1,type='Neyman')# 0.1000000 0.4753221 0.2123389 0.2123389
vector501=dabcd_max_power_test(B=0.1,w=c(1,1,1),p=as.vector((c(85,74,99)+1)/(c(143,138,160)+2) ),nn=c(143,138,160)+2,ub=c(500,500,500)+2,armn=3) # 0.1000000 0.4755912 0.4244088 (0.4244088/2=0.2122044)
all.equal(vector401, c(vector501[1:2],vector501[3]/2,vector501[3]/2),tolerance=1e-3)

#########Five Arms
vector11=dabcd_max_power(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5,BB=0.2, type='Neyman')#0.2 0.2 0.2 0.2 0.2
vector21=dabcd_max_power_test(B=0.2,w=c(1,1,1,1,1),p=as.vector((c(54,67,85,63,70)+1)/(c(100,88,90,94,102)+2) ),nn=c(100,88,90,94,102)+2,ub=c(600,600,600,600,600),armn=5)#0.2 0.2 0.2 0.2 0.2
all.equal(vector11, vector21)

vector31=dabcd_max_power(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5,BB=0.1, type='Neyman')#0.4234212 0.1000000 0.2765788 0.1000000 0.1000000
vector41=dabcd_max_power_test(B=0.1,w=c(1,1,1,1,1),p=(c(54,67,85,63,70)+1)/(c(100,88,90,94,102)+2),nn=c(100,88,90,94,102)+2,ub=c(500,500,500,500,500),armn=5)#0.4234212 0.1000000 0.2765788 0.1000000 0.1000000
all.equal(vector31, vector41)

vector51=dabcd_max_power(NN=c(120,63,87,94,150),Ntotal1=c(200,88,110,108,160),armn=5,BB=0.15, type='Neyman' )#0.3442241 0.1500000 0.1500000 0.1500000 0.2057759
vector61=dabcd_max_power_test(B=0.15,w=c(1,1,1,1,1),p=(c(120,63,87,94,150)+1)/(c(200,88,110,108,160)+2),nn=c(200,88,110,108,160)+2,ub=c(300,300,300,300,300),armn=5)#0.3442241 0.1500000 0.1500000 0.1500000 0.2057759
all.equal(vector51, vector61)

vector71=dabcd_max_power(NN=c(120,120,120,120,30),Ntotal1=c(150,150,150,150,90),armn=5,BB=0.1, type='Neyman' )#0.1150440 0.1150440 0.1150440 0.1150440 0.5398239
p=c(0.7960526,0.7960526, 0.7960526, 0.7960526, 0.3369565)
q=1-p
sqrt(p[5]*q[5])/(sqrt(p[1]*q[1])+sqrt(p[5]*q[5]))#0.5398238
sqrt(p[1]*q[1])/4/(sqrt(p[1]*q[1])+sqrt(p[5]*q[5]))#0.115044
all.equal(vector71, c(0.115044,0.115044,0.115044,0.115044,0.5398238),tolerance=1e-6)

vector81=dabcd_max_power(NN=c(120,120,120,120,30),Ntotal1=c(150,150,150,150,90),armn=5,BB=0.2, type='Neyman')#0.2 0.2 0.2 0.2 0.2
all.equal(vector81, c(0.2,0.2,0.2,0.2,0.2))


