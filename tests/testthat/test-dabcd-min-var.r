set.seed(12345)
vector11=dabcd_min_var(NN=c(25, 25, 25)*c(0.4, 0.28, 0.6),Ntotal1=c(25, 25, 25),armn=3,type='RSIHR',dabcd=1,gamma=2)
all.equal(vector11,c(0.2963430,0.1837979,0.5198591),tolerance=1e-6,scale=1)

vector21=dabcd_min_var(NN=c(82, 79, 88)*c(0.39, 0.43, 0.4),Ntotal1=c(82, 79, 88), armn=3,type='Neyman',dabcd=1,gamma=2)
all.equal(vector21,c(0.3329618,0.3743423,0.2926959),tolerance=1e-6,scale=1)

vector31=dabcd_min_var(NN=c(77,80,83,90,94)* c(0.8, 0.7,0.67,0.78,0.75),Ntotal1=c(77,80,83,90,94),armn=5,type='Neyman',dabcd=1,gamma=2)
all.equal(vector31,c(0.1871102,0.2534742,0.2532257,0.1502448,0.1559451),tolerance=1e-6,scale=1)

vector41=dabcd_min_var(NN=c(77,80,83,90,94)* c(0.8, 0.7,0.67,0.78,0.75),Ntotal1=c(77,80,83,90,94),armn=5,type='RSIHR',dabcd=1,gamma=2)
all.equal(vector41,c(0.2672184,0.2033942,0.1772075,0.1887906,0.1633893),tolerance=1e-6,scale=1)

vector51=dabcd_min_var(NN=c(34,45,65,45,57)*c(0.4, 0.3,0.38,0.3,0.36),Ntotal1=c(34,45,65,45,57),armn=5,type='RSIHR',dabcd=1,gamma=2)
all.equal(vector51,c(0.4324015,0.1637758,0.1088254,0.1637758,0.1312215),tolerance=1e-6,scale=1)

vector61=dabcd_min_var(NN=c(34,45,65,45,57)*c(0.4, 0.3,0.38,0.3,0.36),Ntotal1=c(34,45,65,45,57),armn=5,type='Neyman',dabcd=1,gamma=2)
all.equal(vector61,c(0.3919406,0.1862491,0.1041620,0.1862491,0.1313992),tolerance=1e-6,scale=1)


