#### IFM Simulation wroted by Rostam Abdollahi-Arpanahi ####
############################
#### Input parameters 	####
############################
rm(list=ls())

h2 <- 0.3; PhenVar <- 100; GenVar <- h2*PhenVar;
ErrVar <- PhenVar-GenVar
Mu <- 20

### size of BasePop and no of progenies per generations
nBase <- 100
nPrg <- 100
### Percentage of Sires selected in each generation
sr  = 0.10
### Percentage of Dams selected in each generation
dr  = 0.90

### no. Generations
nGen <- 10

############################################
#										   #
#### 		Body of Program 			####
#										   #
############################################

#### Simulate the founder population
## Store the founder generation data in the BasePop
BasePop <- matrix(0,nrow=nBase, ncol=7)
colnames(BasePop)<- c("ID", "Sire","Dam","Sex","Gen","BV", "Phen")

### Create Phenotypes for founder generation
for(i in 1:nBase) {
	### Breeding values
	BasePop[i,6] <- rnorm(1)*sqrt(GenVar)
	### Residuals
	Residual <- rnorm(1)*sqrt(ErrVar)
	### Phenotypes=mean+ BV+ Error
	BasePop[i,7] <- Mu+ BasePop[i,5]+Residual
	### Gender
	BasePop[i,4] <- sample(1:2,1,TRUE)
	### Animal ID
	BasePop[i,1] <- i
}

## Keep the parent data for generating the progenies
Pop <- BasePop
#### Create a progeny matrix to save progeny records

Progeny <- matrix(0,nrow=nPrg, ncol=7)

### simulate progeny records

for (Gen in 1:nGen) {

### Extract sires from the data and sort based on True Breeding Values
Sire_data <- Pop[Pop[,4]==1,]
Sire_data <- Sire_data[order(Sire_data [,6],decreasing = TRUE),]

### Extract dams from the data and sort based on True Breeding Values
Dam_data <-  Pop[Pop[,4]==2,]
Dam_data <- Dam_data[order(Dam_data[,6],decreasing = TRUE),]

## no. of sires for matings
no.sire <- sr*nrow(Sire_data)
## select sire IDs
sires.list <- Sire_data[1:no.sire,1]

## no. of dams for mating
no.dam <- dr*nrow(Dam_data)
# select dam IDs
dams.list <-Dam_data[1:no.dam,1]
 
### Loop for generating progenies
	for(prg in 1:nPrg){
		
		### sample sire and dam as the perent of progeny
		sire <- sample(sires.list,1)
		dam <- sample(dams.list,1)
		
		### Generating the breeding values for progeny
		Progeny[prg,6] <- 0.5*(Pop[sire,6]+Pop[dam,6])+(sqrt(0.5*GenVar)*rnorm(1))
		### Phenotype for progenies
		Progeny[prg,7] <- Mu+ Progeny[(prg),6]+ rnorm(1)*sqrt(ErrVar)
		
		## Asign sire, dam, Gender and generation number to the progeny
		### sex
		Progeny[prg,4] <-  sample(1:2,1,TRUE)
		### sir identification
		Progeny[prg,2] <-  sire
		### dam identification
		Progeny[prg,3] <-  dam
		### Generation number
		Progeny[prg,5] <-  Gen
		
		### Animal IDs
		Progeny[prg,1] <-  prg+ nrow(Pop)	
	
	
	}
	
	### Append progenies to the old data file (Pop)
	Pop <- rbind(Pop,Progeny)
}


### Save data file

write.table(Pop , file = "data_ST.txt", append = FALSE, quote = FALSE, sep = " ",row.names = FALSE)

### Save pedigree file

Ped <- Pop[,1:3]
colnames(Ped) <- c("Progeny", "sire", "dam")
write.table(Ped , file = "ped_ST.txt", append = FALSE, quote = FALSE, sep = " ",row.names = FALSE)




