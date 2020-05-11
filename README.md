# Single and Multiple-trait Infinitesimal Model Simulation

**Author**: *Rostam Abdollahi-Arpanahi*

**Date**: May 10th, 2020

---

- Many animal and plant breeding students are looking for a simple R code to learn the stochastic simulation. Two R-scripts for single- and multiple-trait infinitesimal simulation have been provided here.  
- In the Multiple trait simulation, seven traits are simulated simultaneously and the selection of the best animals is based on the selection index. Hence, the reader should be familiar with the selection index theory to be able to digest the code. The Genetic and Phenotypic (co)varaince matrices between these seven traits are given in the beginning of script.
- These two scripts are excellent resources for teaching in the class and also for learning the base of infinitesimal simulation. In case, you are interested in finite loci model simulation, I would suggest you use QMSim program written by Dr. Mehdi Sargolzaei. Download it from [here](http://animalbiosciences.uoguelph.ca/~msargol/qmsim/).

---

- Run the <u>*Single trait*</u> infinitesimal model simulation

```
> git clone https://github.com/Rostamabd/Single-and-Multi-trait-IFM-Simulation.git
> module load R
> R
> source("Single-Trait-Simulation.R")
```

- output files

  - data_ST.txt: This file contains "Animal_ID, Sire, Dam, Gender, Generation, BV, Phenotype".

  - ped_ST.txt: This file contains "Animal_ID, Sire, Dam"

    

- Run the <u>*Multi-trait*</u> infinitesimal model simulation 

```
> git clone https://github.com/Rostamabd/Single-and-Multi-trait-IFM-Simulation.git
> module load R
> R
> source("Multi_trait_Simulation.R")
```

- output files
  - data_MT.txt: This file contains "ID, Sire, Dam, Gender, BV1, BV2, BV3, BV4, BV5, BV6, BV7, aggregate genotype(H) and Index".
  - ped_MT.txt: This file contains "Animal_ID, Sire and Dam, Generation"

## Contact Information

Please send your comments and suggestions to rostam7474 at gmail dot com