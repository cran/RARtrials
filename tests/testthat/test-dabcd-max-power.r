set.seed(12345)
#####Three-Armed DABCE Targeting Neyman Allocation

vector11=dabcd_max_power(NN=c(17,9,19),Ntotal1=c(129,155,123),armn=3, BB=0.25, type='Neyman',dabcd=0)
all.equal(vector11,c(0.2500000,0.3688891,0.3811109),tolerance=1e-6,scale=1)

vector21=dabcd_max_power(NN=c(1,1,2),Ntotal1=c(3,5,4),armn=3, BB=0.3, type='Neyman',dabcd=0)
all.equal(vector21,c(0.3000000,0.3488257,0.3511743),tolerance=1e-6,scale=1)

vector31=dabcd_max_power(NN=c(3, 6, 8),Ntotal1=c(34, 36, 34),armn=3, BB=0.2,type='Neyman',dabcd=0)
all.equal(vector31, c(0.3555952,0.2000000,0.4444048),tolerance=1e-6,scale=1)

vector41=dabcd_max_power(NN=c(30, 20, 18),Ntotal1=c(90, 100, 105),armn=3, BB=0.1, type='Neyman',dabcd=0)
all.equal(vector41, c(0.5274177,0.1000000,0.3725823),tolerance=1e-6,scale=1)

vector51=dabcd_max_power(NN=c(35,66,50),Ntotal1=c(100,88,90),armn=3,BB=0.1, type='Neyman',dabcd=0)
all.equal(vector51, c(0.4698286,0.4301714,0.1000000),tolerance=1e-6,scale=1)

vector61=dabcd_max_power(NN=c(35,66,66),Ntotal1=c(100,90,90),armn=3,BB=0.1, type='Neyman',dabcd=0) 
all.equal(vector61, c(0.5178969,0.2410515,0.2410515),tolerance=1e-6,scale=1)

vector71=dabcd_max_power(NN=c(35,66,66),Ntotal1=c(100,90,90),armn=3,BB=0.25, type='Neyman',dabcd=0)
all.equal(vector71, c(0.50,0.25,0.25),tolerance=1e-6,scale=1)

vector81=dabcd_max_power(NN=c(120,120,30),Ntotal1=c(150,150,90),armn=3,BB=0.1, type='Neyman',dabcd=0) 
all.equal(vector81, c(0.2300881,0.2300881,0.5398239),tolerance=1e-6,scale=1)

vector91=dabcd_max_power(NN=c(120,120,30)+1,Ntotal1=c(150,150,90)+2,armn=3,BB=0.25, type='Neyman',dabcd=0)
all.equal(vector91, c(0.25,0.25,0.50),tolerance=1e-6,scale=1)

#####Three-Armed DABCE Targeting RSIHR Allocation
vector101=dabcd_max_power(NN=c(40,52,66),Ntotal1=c(102,130,130),armn=3,BB=0.1, type='RSIHR',dabcd=0 )
all.equal(vector101, c(0.375952,0.100000,0.524048),tolerance=1e-6,scale=1)

vector102=dabcd_max_power(NN=c(80,38,55),Ntotal1=c(110,88,96),armn=3,BB=0.2, type='RSIHR',dabcd=0  )
all.equal(vector102, c(0.4583773,0.3416227,0.2000000),tolerance=1e-6,scale=1)

vector103=dabcd_max_power(NN=c(80,55,55),Ntotal1=c(110,88,88),armn=3,BB=0.15, type='RSIHR' ,dabcd=0 )
all.equal(vector103, c(0.5187922,0.2406039,0.2406039),tolerance=1e-6,scale=1)

vector104=dabcd_max_power(NN=c(80,55,55),Ntotal1=c(110,88,88),armn=3,BB=0.25, type='RSIHR')
all.equal(vector104, c(0.50,0.25,0.25),tolerance=1e-6,scale=1)

#########Four-Armed DABCE Targeting Neyman Allocation
vector201=dabcd_max_power(NN=c(54,67,85,63),Ntotal1=c(100,88,90,94),armn=4,BB=0.2, type='Neyman')
all.equal(vector201, c(0.3466227,0.2000000,0.2533773,0.2000000),tolerance=1e-6,scale=1)

vector202=dabcd_max_power(NN=c(120,80,111,93),Ntotal1=c(154,128,140,160),armn=4,BB=0.1, type='Neyman')
all.equal(vector202, c(0.1000000,0.1000000,0.3463452,0.4536548),tolerance=1e-6,scale=1)

vector203=dabcd_max_power(NN=c(80,80,80,93),Ntotal1=c(154,154,154,160),armn=4,BB=0.1, type='Neyman')
all.equal(vector203, c(0.1676923,0.1676923,0.1676923,0.4969231),tolerance=1e-6,scale=1)

vector204=dabcd_max_power(NN=c(80,80,80,93),Ntotal1=c(154,154,154,160),armn=4,BB=0.2, type='Neyman')
all.equal(vector204, c(0.2,0.2,0.2,0.4),tolerance=1e-6,scale=1)

vector205=dabcd_max_power(NN=c(54,54,85,63)+1,Ntotal1=c(88,88,90,94)+2,armn=4,BB=0.2, type='Neyman')
all.equal(vector205, c(0.2326693,0.2326693,0.3346614,0.2000000),tolerance=1e-6,scale=1)

vector206=dabcd_max_power(NN=c(154,154,185,163),Ntotal1=c(188,188,190,194),armn=4,BB=0.2,type='Neyman')
all.equal(vector206, c(0.2487669,0.2487669,0.3024662,0.2000000),tolerance=1e-6,scale=1)

vector207=dabcd_max_power(NN=c(120,120,120,30),Ntotal1=c(150,150,150,90),armn=4,BB=0.1, type='Neyman' )
all.equal(vector207, c(0.1533920,0.1533920,0.1533920,0.5398239),tolerance=1e-6,scale=1)

vector208=dabcd_max_power(NN=c(120,120,120,30),Ntotal1=c(150,150,150,90),armn=4,BB=0.25, type='Neyman' )
all.equal(vector208, c(0.25,0.25,0.25,0.25),tolerance=1e-6,scale=1)

vector209=dabcd_max_power(NN=c(85,74,99,99)+1,Ntotal1=c(143,138,160,160)+2,armn=4,BB=0.1,type='Neyman')
all.equal(vector209, c(0.1000000,0.4753221,0.2123389,0.2123389),tolerance=1e-6,scale=1)


#########Five-Armed DABCE Targeting Neyman Allocation
vector301=dabcd_max_power(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5,BB=0.2, type='Neyman')
all.equal(vector301, c(0.2,0.2,0.2,0.2,0.2),tolerance=1e-6,scale=1)

vector302=dabcd_max_power(NN=c(54,67,85,63,70),Ntotal1=c(100,88,90,94,102),armn=5,BB=0.1, type='Neyman')
all.equal(vector302, c(0.4234212,0.1000000,0.2765788,0.1000000,0.1000000),tolerance=1e-6,scale=1)

vector303=dabcd_max_power(NN=c(120,63,87,94,150),Ntotal1=c(200,88,110,108,160),armn=5,BB=0.15, type='Neyman' )
all.equal(vector303, c(0.3442241,0.1500000,0.1500000,0.1500000,0.2057759),tolerance=1e-6,scale=1)

vector304=dabcd_max_power(NN=c(120,120,120,120,30),Ntotal1=c(150,150,150,150,90),armn=5,BB=0.1, type='Neyman' )
all.equal(vector304, c(0.1150440,0.1150440,0.1150440,0.1150440,0.5398239),tolerance=1e-6,scale=1)

vector305=dabcd_max_power(NN=c(120,120,120,120,30),Ntotal1=c(150,150,150,150,90),armn=5,BB=0.2, type='Neyman')
all.equal(vector305, c(0.2,0.2,0.2,0.2,0.2),tolerance=1e-6,scale=1)



