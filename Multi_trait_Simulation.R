#### Multiple trait Infinitesimal Simulation ####
#################################
#### Input parameters 		  ###
#################################

rm(list=ls())
library(Matrix)

## all studied traits
all_traits=7

### Traits included in the Selection Index (I)
ntrt <- 4
TrtNum <- c(1,2,3,4)

### Genetic variance_Covariance
cov_gen <- matrix(c(147.35,	39.63991319,43.59509248,20.27031019,	0.11,0.312452219,0.112162335,
39.63991319,31.7,17.62813932,5.062553703,0.029727854,0.119458649,0.054773211,
43.59509248,17.62813932,21.2,7.688706003,0.076595979,0.167035503,0.040518243,
20.27031019,5.062553703,7.688706003,6.6,0.018086087,0.170848607,0.028259512,
0.11,0.029727854,0.076595979,0.018086087,0.001936,0.0031,0.0016,
0.312452219,0.119458649,0.167035503,0.170848607,0.0031,	0.0098,	0.00544,
0.112162335,0.054773211,0.040518243,0.028259512,0.0016,	0.00544,0.0084
),nrow=7,byrow=T)

## Convert the Genetic (co)variance matrix to positive definite
gg <- nearPD(cov_gen)

## New Positive definite Genetic (co)variance matrix
cov_gen <- gg$mat

## Genetic (co)variance matrix between traits in Index and traits in Breeding Objective
cov_gen2  <- matrix(cov_gen [TrtNum,c(1,5,6,7)],nrow=ntrt, byrow=F)

### Phenotypic Varinace_Covariance
cov_phen <- matrix(c(421,52.64529609,	45.04800466,	30.10388627,	2,	2.2,	3.53,
52.64529609,	107.02,	17.62813932,	8.307395801,	0.18103,	0.15517,	0.6362,
45.04800466,	10.628627,	118.02,	8.430622753,	0.76046,	0.48887,	0.890823,
30.10388627,	8.307395801,	8.430622753,	34.64,	0.618,	0.53,	0.120654,
2,	0.18103,	0.76046,	0.618,	0.1225,	0.063,	0.001613,
2.2,	0.15517,	0.48887,	0.53,	0.063,	0.09,	0.003629,
3.53,	0.6362,	0.890823,	0.120654,	0.001613,	0.003629,	0.1681),nrow=7,byrow=T)

## Phenotypic (co)variance matrix between traits in Index
cov_phen2  <- matrix(cov_phen [TrtNum,TrtNum],nrow=ntrt, byrow=F)

### Environmental (co)variance matrix
Cov_Env <- cov_phen-cov_gen

### Convert environmental (co)variance matrix to nearest positive definite
ee <- nearPD(Cov_Env)

## Positive definit environmental (co)variance matrix
Cov_Env <- ee$mat

### Trait means
Mu <- c(185,63.5,102.5,104.8,0.96,0.89,0.96)

#### Genetic(co)varinaces between traits in Selection Index
gen_Ind_cov <- cov_gen [TrtNum,TrtNum]


#### Economic values of traits included in the Breeding Objective 
Ecnm_val <- c(37735,39722,96405,39722)

### Partial regression coefficients of selection Index (I=b1*x1+b2*x2+b3*x3+b4*x4)
b <- solve(cov_phen2)%*% (cov_gen2%*%Ecnm_val)

### Size of Base Population
nPop <- 100

### number of progenies per generation
nPrg <- 100
### Number of Generations
nGen <- 10

### Percentage of Sires selected in each generation
sr  = 0.10
### Percentage of Dams selected in each generation
dr  = 0.90
### Cholesky decomposition of (co)variance  martices 

Lg <- t(chol(cov_gen))

Le <- t(chol(Cov_Env))
#############################################
#### 		2. Body of Program 		      ###
#############################################


#### Founder population			#############

## Create a jar for storing the base population

BasePop <- matrix(0,nrow=nPop, ncol=(6+all_traits))
ff =paste("BV",1:all_traits,sep="")
colnames(BasePop)<- c("ID", "Sire","Dam","Sex",ff,"H", "Index")


### Create Phenotypes for Base Population
for(i in 1:nPop) {
	# generate the breeding values(BV) for seven traits
	Bv <- Lg%*%matrix(rnorm(all_traits),nrow=all_traits,ncol=1)
	# generate the environmental effects for seven traits
	Ev <- Le%*%matrix(rnorm(all_traits),nrow=all_traits,ncol=1)
	# Phenotype is equal to mean+BV+EV
	pheno <- Mu+Bv+Ev
	# Total aggeragte genotypic value (Agg_G)
	Agg_G <- t(Ecnm_val)%*%Bv[c(1,5,6,7)]
	#The value of Index (Index) per each individual
	Index <- t(b)%*%pheno[TrtNum]
	# Assign variables to the Base population matrix
	BasePop[i,5:(5+all_traits-1)] <- Bv[1:all_traits]
	BasePop[i,(5+all_traits)] <- Agg_G
	BasePop[i,(6+all_traits)] <- Index
	
}

# The first four columnas are ID, Sire, Dam and gender
BasePop[,1] <- 1:nPop
BasePop[,2]	<- 0
BasePop[,3] <- 0
BasePop[,4] <- sample(1:2,nPop,TRUE)

### Form the pedigree of base population
Ped <- cbind(BasePop[,1:3],0)
colnames(Ped) <- c("Progeny","Sire","Dam","Gen")

### Store all the base and progeny into a new jar called Pop
Pop <- BasePop


## Herd is a temporary matrix keeping the parent generation
Herd <- BasePop


#### Non-founder generations (non-overlapping generations)
Counter <- nPop
for(Gen in 1:nGen){

	count <- table(Herd[, 4])

	count <-as.vector(count)

	## Extract list of Males and Females
	List_Sire <- as.matrix(subset(Herd, Herd[, 4] == 1))
	List_Dam <- as.matrix(subset(Herd , Herd[, 4] == 2))

	### Sort Sires and Dams list based on the Index
		
	Sire.list <- List_Sire[order(List_Sire[,(6+all_traits)],decreasing = TRUE),]
		
	Dam.list <- List_Dam[order(List_Dam[,(6+all_traits)],decreasing = TRUE),]

	# Generate the breeding values for progenies	
	Progeny <- matrix(NA,nPrg,(6+all_traits))

		for(i in 1:nPrg){
			Bv <- as.numeric()
			Ev <- as.numeric()
			
			## Here, select randomly dams
			no.dam <- dr*nrow(Dam.list)
			dam.select <- sample(Dam.list[1:no.dam,1],1)

			## Here, select 10 percenatge of top sires
			no.sire <- sr*nrow(Sire.list)
			sire.select <- sample(Sire.list[1:no.sire,1],1)
			## Generate the breeding values
			for(t in 1:all_traits){
					Bv[t] <-0.5*(as.vector(na.omit(Sire.list[Sire.list[, 1]==sire.select,(4+t)]))+Herd[Herd[, 1]==dam.select,(4+t)])
			
				}
			## Mendelian sampling
			ms <- sqrt(0.5)*(Lg%*%matrix(rnorm(all_traits),nrow=all_traits,ncol=1))
			# Total Breeding value is equal to BV plus mendelian sampling
			Bv <- Bv+ ms
			# Add random environmental effect to the breeding values
			Ev <- Le%*%matrix(rnorm(all_traits),nrow=all_traits,ncol=1)
			#generate Phenotypic value for progenies
			pheno <- Mu[1:all_traits]+Bv+Ev
			# Aggregate genotypic value for progenies (it is recommanded to recalculate economic values in each generation)
			Agg_G <- t(Ecnm_val)%*%Bv[c(1,5,6,7)]
			# Total Index per each progeny (it is recommanded to recalculate b values in each generation)
			Index <- t(b)%*%pheno[TrtNum]
			
			# Asing BV, Index, H and pedigree information to the progenies
			Progeny[i,5:(5+all_traits-1)] <- Bv[1:all_traits]
			Progeny[i,(5+all_traits)] <- Agg_G
			Progeny[i,(6+all_traits)] <- Index
			
			### Define the new animal IDs
			Counter= Counter+1
			Progeny[i,1] <- Counter
			# Sire ID
			Progeny[i,2]	<- sire.select
			#Dam ID
			Progeny[i,3] <- dam.select
			#Gender
			Progeny[i,4] <- sample(1:2,1,TRUE)
		}
		
	Herd <- Progeny
	## store all the data and pedigree together
	Pop <- rbind(Pop,Progeny)
	Ped <- rbind(Ped,cbind(Progeny[,1:3],Gen))
}

## Add generation number to the data, in order to check the genetic and phenotypic trends
Pop <- round(cbind(Pop,Ped[,4]),digits = 3)

## Save the data and pedigree
write.table(Pop ,file="data_MT.txt",quote = F,row.names=F)

write.table(Ped,file="ped_MT.txt",quote = F,row.names=F)
	
	