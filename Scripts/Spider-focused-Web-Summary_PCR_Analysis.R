## USES ALL SPIDERS (i.e. does not exclude hungry spiders)
## Change working directory, dataset required, and settings below, and then run the entire script.
## The number of bootstrap replications must be a multiple of 100


## Requires the library sqldf, which requires the following:
	## gsubfn, proto, RSQLite, DBI, tcltk2, chron

## Requires the library bipartite, which requires the following:
	## vegan, permute, lattice
	## sna, network, statnet.common, igraph,
	## fields, spam, maps, mapproj

## Requires the library dostats (most recent version built under R version 3.2.1)

## Requires the library matrixStats (most recent version built under R version 3.2.1)

## Requires the library exact2x2 (most recent version built under R version 3.2.1), which requires the following:
	## exactci
	## ssanv

## Program Returns the following information:
	## Percentile Bootstrap Confidence Intervals
		## Undirected Connectance with Cannibalism
		## Undirected Connectance
		## Directed Connectance
		## Interaction Strength
		## IGP Undirected Connectance Predators
		## IGP Directed Connectance Predators
		## IGP Directed Connectance Prey
		## IGP Interaction Strength Predators
		## IGP Interaction Strength Prey
		## Interaction Evenness; logbase = e; intereven = sum
		## Interaction Evenness; logbase = 2; intereven = sum
		## Interaction Evenness; logbase = e; intereven = prod
		## Interaction Evenness; logbase = 2; intereven = prod
		## Specialization dprime; Pred = col, Prey = row
			## na.rm If TRUE, NAs are excluded first, otherwise not
			## drop If TRUE, singleton dimensions in the result are dropped, otherwise not
			## Settings: na.rm = FALSE, drop = FALSE
		## Specialization dprime transpose; Prey = col, Pred = row
			## na.rm If TRUE, NAs are excluded first, otherwise not
			## drop If TRUE, singleton dimensions in the result are dropped, otherwise not
			## Settings: na.rm = FALSE, drop = FALSE
		## Specialization h2; Pred = col, Prey = row
		## Specialization h2 uncorrelated; Pred = col, Prey = row
		## Specialization h2 transpose; Prey = col, Pred = row
		## Specialization h2 uncorrelated transpose; Prey = col, Pred = row
		## Compartmentalization from Bipartite
		## Compartmentalization Pimm & Lawton; All Spiders
		## Compartmentalization Pimm & Lawton; Predators Only


	## Median Values
		## Undirected Connectance with Cannibalism
		## Undirected Connectance
		## Directed Connectance
		## Interaction Strength
		## IGP Undirected Connectance Predators
		## IGP Directed Connectance Predators
		## IGP Directed Connectance Prey
		## IGP Interaction Strength Predators
		## IGP Interaction Strength Prey
		## Interaction Evenness; logbase = e; intereven = sum
		## Interaction Evenness; logbase = 2; intereven = sum
		## Interaction Evenness; logbase = e; intereven = prod
		## Interaction Evenness; logbase = 2; intereven = prod
		## Specialization dprime; Pred = col, Prey = row
			## na.rm = If TRUE, NAs are excluded first, otherwise not
			## Settings: na.rm = FALSE
		## Specialization dprime transpose; Prey = col, Pred = row
			## na.rm = If TRUE, NAs are excluded first, otherwise not
			## Settings: na.rm = FALSE
		## Specialization h2; Pred = col, Prey = row
		## Specialization h2 uncorrelated; Pred = col, Prey = row
		## Specialization h2 transpose; Prey = col, Pred = row
		## Specialization h2 uncorrelated transpose; Prey = col, Pred = row
		## Compartmentalization from Bipartite
		## Compartmentalization Pimm & Lawton; All Spiders
		## Compartmentalization Pimm & Lawton; Predators Only

	## Where data_type = ALL results
		## Contingency Table Results Quant Data Median Odds Ratios
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Quant Data Median P-Values
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Indicator Data Median Odds Ratios
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Indicator Data Median P-Values
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Quant Data - Confidence Interval on the Odds Ratios
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Quant Data - Confidence Interval on the P-Values
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Indicator Data - Confidence Interval on the Odds Ratios
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test
		## Contingency Table Results Indicator Data - Confidence INterval on the P-Values
			## Central Fisher Exact Test
			## Fisher's Exact Test with Minimum Likelihood
			## Blaker's Exact Test


	## Original Dataset
		## Undirected Connectance Value with Cannibalism
		## Undirected Connectance Value
		## IGP Predators Undirected Connectance Value

		## Contingency Table for All Spiders Only
			## Original Dataset Contingency Table Indicator Results
				## Central Fisher Exact Test
				## Fisher's Exact Test with Minimum Likelihood
				## Blaker's Exact Test
			## Original Dataset Contingency Table Quant Results
				## Central Fisher Exact Test
				## Fisher's Exact Test with Minimum Likelihood
				## Blaker's Exact Test



## Notes:
	## 1. Runs EXTRA Results for Original Dataset only.  Must refer to other R programs for original dataset.
	## 2. Contingency Table only generates when run for All Spiders.  Cursorial only and Web only will not generate these results.
	## 3. For ALL Spiders only, Confidence Intervals for each bootstrap run are
		## stored in results.xxxx_xxxx_confint.  Results are NOT automatically generated.
		## To generate Results, uncomment print section.

#load packages
#library
library(dplyr)
library(tidyverse)

#=====import data======
#relative pathname
pcr <- file.path(".", "Data", "Summary_PCR_RawData.csv")
print(pcr)

#import data
dataset<- read_csv(pcr)

## Settings
alpha <- 0.05				## 95% Confidence Interval
B <- 5000					## Number of Bootstrap Replications
n <- 299					## Number Spiders to Select per Family
numberpredatoritems <- 11		## For IntraGuildPredation Number Columns of Predators
numberpreyitems <- 11			## For IntraGuildPredation Number Columns of Prey
data_type = "ALL"				## Data type: ALL = All Spiders, CURSORIAL = Cursorial Spiders,
							## WEB = Web Building Spiders






######################################### DO NOT MODIFY BELOW THIS LINE
## Extra Libraries Needed
	library(sqldf)
	library(bipartite)
	library(dostats)
	library(matrixStats)
	library(exact2x2)
	library(tcltk2)




## Blank Datasets Required (do not adjust)
spiders.selected.agel <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.cori <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.dict <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.gnap <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.hahn <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.liny <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.lyco <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.pisa <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.salt <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.ther <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
spiders.selected.thom <- matrix(, nrow = B, ncol = n, byrow = TRUE)				## Storing Spiders Selected during each run
results.Connectance_wcannibalism <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Undirected Connectance with Cannibalism Results Vector for All Species
results.Connectance_undirected <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Undirected Connectance Results Vector for All Species
results.Connectance_directed <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Directed Connectance Results Vector
results.Interaction_Strength <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Interaction Strength Results Vector
results.IGP_Conn_Pred_undirected <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Undirected Connectance Results Vector for IGP Connectance Predators
results.IGP_Conn_Pred_directed <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Directed Connectance Results Vector for IGP Connectance Predators
results.IGP_Conn_Prey_directed <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Directed Connectance Results Vector for IGP Connectance Prey
results.IGP_IS_Pred <- matrix(, nrow = 1, ncol = B, byrow = TRUE)					## Interaction Strength Results Vector for IGP IS Predators
results.IGP_IS_Prey <- matrix(, nrow = 1, ncol = B, byrow = TRUE)					## Interaction Strength Results Vector for IGP IS Prey
results.Compart_Bipartite <- matrix(, nrow = 1, ncol = B, byrow = TRUE)				## Number of Compartments using Bipartite
results.Interaction_Evenness_e_sum <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Interaction Evenness Results logbase=e intereven=sum
results.Interaction_Evenness_2_sum <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Interaction Evenness Results logbase=2 intereven=sum
results.Interaction_Evenness_e_prod <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Interaction Evenness Results logbase=e intereven=prod
results.Interaction_Evenness_2_prod <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Interaction Evenness Results logbase=2 intereven=prod
results.Specialization_dprime <- matrix(data = NA, 
	nrow = (numberpredatoritems + numberpreyitems),
	ncol = B, byrow = TRUE)											## Specialization dprime Results Matrix
	rownames(results.Specialization_dprime) <- colnames(dataset[3:24])
results.Specialization_h2 <- matrix(, nrow = 1, ncol = B, byrow = TRUE)				## Specialization h2 Results Vector
results.Specialization_h2_uncorr <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Specialization h2 uncorrelated Results Vector
results.Specialization_pred_h2 <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Specialization predators h2 Results Vector
results.Specialization_pred_h2_uncorr <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Specialization predators h2 uncorrelated Results Vector
results.Compartmentalization_pimm_all <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Compartmentalization Pimm & Lawton Results Vector All Spiders
results.Compartmentalization_pimm_pred <- matrix(, nrow = 1, ncol = B, byrow = TRUE)	## Compartmentalization Pimm & Lawton Results Vector Predators / Predators
if(data_type == "ALL"){
	results.quant_fisher_central_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)## Quant Data; Central Fisher's Exact Test Odds Ratio Vector
	results.quant_fisher_central_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)	## Quant Data; Central Fisher's Exact Test P-value Vector
	results.quant_fisher_central_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Quant Data; Central Fisher's Exact Test Confidence Intervals Matrix

	results.quant_fisher_minlike_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)## Quant Data; Fisher Exact Test with Minimum Likelihood Odds Ratio Vector
	results.quant_fisher_minlike_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)	## Quant Data; Fisher Exact Test with Minimum Likelihood P-value Vector
	results.quant_fisher_minlike_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Quant Data; Fisher Exact Test with Minimum Likelihood Confidence Intervals Matrix

	results.quant_blaker_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Quant Data; Blaker Exact Test Odds Ratio Vector
	results.quant_blaker_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Quant Data; Blaker Exact Test P-value Vector
	results.quant_blaker_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Quant Data; Blaker Exact Test Confidence Intervals Matrix

	results.ind_fisher_central_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)	## Indicator Data; Central Fisher's Exact Test Odds Ratio Vector
	results.ind_fisher_central_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Indicator Data; Central Fisher's Exact Test P-value Vector
	results.ind_fisher_central_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Indicator Data; Central Fisher's Exact Test Confidence Intervals Matrix

	results.ind_fisher_minlike_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)	## Indicator Data; Fisher Exact Test with Minimum Likelihood Odds Ratio Vector
	results.ind_fisher_minlike_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Indicator Data; Fisher Exact Test with Minimum Likelihood P-value Vector
	results.ind_fisher_minlike_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Indicator Data; Fisher Exact Test with Minimum Likelihood Confidence Intervals Matrix

	results.ind_blaker_oddsratio <- matrix(, nrow = 1, ncol = B, byrow = TRUE)		## Indicator Data; Blaker Exact Test Odds Ratio Vector
	results.ind_blaker_pval <- matrix(, nrow = 1, ncol = B, byrow = TRUE)			## Indicator Data; Blaker Exact Test P-value Vector
	results.ind_blaker_confint <- matrix(data = NA,
		nrow = B, ncol = 2, byrow = TRUE)								## Indicator Data; Blaker Exact Test Confidence Intervals Matrix
	}


vec0 <- rep(0,times=22)												## Row with All Zero's
Link.yes <- 1													## Used when going from quantitative dataset to indicator dataset
Link.no <- 0													## Also used later with original dataset calculations
pb_boot <- tkProgressBar(title = "bootstrap progress bar", min = 0, max = B, width = 300)	## Progress Bar to for bootstrap



## Generate Restricted Datasets
		## Count Number of Spiders per Family
		for (i in colnames(dataset)) {
			sql1 <- fn$identity("
				Select Family_Short,
					count(ID) as Total_Spiders
				From dataset
				Group by Family_Short ")
				Num_Per_Fam <- sqldf(c(sql1))
			}

		## Generate a List of the Remaining Families and Counts
		for (i in colnames(dataset)) {
			sql1 <- fn$identity("
				Select dataset.Family_Short,
					count(dataset.ID) as Total_Spiders
				From dataset, Num_Per_Fam
				Where dataset.Family_Short = Num_Per_Fam.Family_Short
					And Num_Per_Fam.Total_Spiders >= '$n'
				Group by dataset.Family_Short ")
				Num_Per_Fam_Restricted <- sqldf(c(sql1))
			}

		## Split Each Family into a Different Table	## Note: Could have used subset() function
		for (i in colnames(dataset)) {
			sql1 <- fn$identity("Select * From dataset Where Family_Short = 'Agel' ")
			sql2 <- fn$identity("Select * From dataset Where Family_Short = 'Cori' ")
			sql3 <- fn$identity("Select * From dataset Where Family_Short = 'Dict' ")
			sql4 <- fn$identity("Select * From dataset Where Family_Short = 'Gnap' ")
			sql5 <- fn$identity("Select * From dataset Where Family_Short = 'Hahn' ")
			sql6 <- fn$identity("Select * From dataset Where Family_Short = 'Liny' ")
			sql7 <- fn$identity("Select * From dataset Where Family_Short = 'Lyco' ")
			sql8 <- fn$identity("Select * From dataset Where Family_Short = 'Pisa' ")
			sql9 <- fn$identity("Select * From dataset Where Family_Short = 'Salt' ")
			sql10 <- fn$identity("Select * From dataset Where Family_Short = 'Ther' ")
			sql11 <- fn$identity("Select * From dataset Where Family_Short = 'Thom' ")
			dataset_Agel <- sqldf(c(sql1))
			dataset_Cori <- sqldf(c(sql2))
			dataset_Dict <- sqldf(c(sql3))
			dataset_Gnap <- sqldf(c(sql4))
			dataset_Hahn <- sqldf(c(sql5))
			dataset_Liny <- sqldf(c(sql6))
			dataset_Lyco <- sqldf(c(sql7))
			dataset_Pisa <- sqldf(c(sql8))
			dataset_Salt <- sqldf(c(sql9))
			dataset_Ther <- sqldf(c(sql10))
			dataset_Thom <- sqldf(c(sql11))
			}

## Set Seed
set.seed(16)

## Select New Spiders for Datasets
for(j in 1:B) {
	## Sample Spiders (each row will be a separate bootstrap run)
		## Generate Separate Datasets for Each Family
		if(Num_Per_Fam_Restricted$Family_Short %contains% "Agel") {
			iz_Agel <- sample(1:nrow(dataset_Agel), n, replace = TRUE)
			spiders.selected.agel[j,] <- t(iz_Agel)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Cori") {
			iz_Cori <- sample(1:nrow(dataset_Cori), n, replace = TRUE)
			spiders.selected.cori[j,] <- t(iz_Cori)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Dict") {
			iz_Dict <- sample(1:nrow(dataset_Dict), n, replace = TRUE)
			spiders.selected.dict[j,] <- t(iz_Dict)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Gnap") {
			iz_Gnap <- sample(1:nrow(dataset_Gnap), n, replace = TRUE)
			spiders.selected.gnap[j,] <- t(iz_Gnap)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Hahn") {
			iz_Hahn <- sample(1:nrow(dataset_Hahn), n, replace = TRUE)
			spiders.selected.hahn[j,] <- t(iz_Hahn)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Liny") {
			iz_Liny <- sample(1:nrow(dataset_Liny), n, replace = TRUE)
			spiders.selected.liny[j,] <- t(iz_Liny)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Lyco") {
			iz_Lyco <- sample(1:nrow(dataset_Lyco), n, replace = TRUE)
			spiders.selected.lyco[j,] <- t(iz_Lyco)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Pisa") {
			iz_Pisa <- sample(1:nrow(dataset_Pisa), n, replace = TRUE)
			spiders.selected.pisa[j,] <- t(iz_Pisa)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Salt") {
			iz_Salt <- sample(1:nrow(dataset_Salt), n, replace = TRUE)
			spiders.selected.salt[j,] <- t(iz_Salt)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Ther") {
			iz_Ther <- sample(1:nrow(dataset_Ther), n, replace = TRUE)
			spiders.selected.ther[j,] <- t(iz_Ther)
			}

		if(Num_Per_Fam_Restricted$Family_Short %contains% "Thom") {
			iz_Thom <- sample(1:nrow(dataset_Thom), n, replace = TRUE)
			spiders.selected.thom[j,] <- t(iz_Thom)
			}
	} ## END SELECTION OF SPIDERS FOR BOOT DATASETS

## Cleanup From Above
	if(exists("iz_Agel") == TRUE){rm(iz_Agel)}
	if(exists("iz_Cori") == TRUE){rm(iz_Cori)}
	if(exists("iz_Dict") == TRUE){rm(iz_Dict)}
	if(exists("iz_Gnap") == TRUE){rm(iz_Gnap)}
	if(exists("iz_Hahn") == TRUE){rm(iz_Hahn)}
	if(exists("iz_Liny") == TRUE){rm(iz_Liny)}
	if(exists("iz_Lyco") == TRUE){rm(iz_Lyco)}
	if(exists("iz_Pisa") == TRUE){rm(iz_Pisa)}
	if(exists("iz_Salt") == TRUE){rm(iz_Salt)}
	if(exists("iz_Ther") == TRUE){rm(iz_Ther)}
	if(exists("iz_Thom") == TRUE){rm(iz_Thom)}

## Generate Results for Each Bootstrap Run
	for(j in 1:B){
		## Create Bootstrap Table
			boot <- data.frame()												## Bootstrap Dataset

		## Put Sampled Spiders in a New Dataset
			## Generate Separate Datasets for Each Family
			if(Num_Per_Fam_Restricted$Family_Short %contains% "Agel") {
				boot_Agel <- dataset_Agel[spiders.selected.agel[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Cori") {
				boot_Cori <- dataset_Cori[spiders.selected.cori[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Dict") {
				boot_Dict <- dataset_Dict[spiders.selected.dict[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Gnap") {
				boot_Gnap <- dataset_Gnap[spiders.selected.gnap[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Hahn") {
				boot_Hahn <- dataset_Hahn[spiders.selected.hahn[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Liny") {
				boot_Liny <- dataset_Liny[spiders.selected.liny[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Lyco") {
				boot_Lyco <- dataset_Lyco[spiders.selected.lyco[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Pisa") {
				boot_Pisa <- dataset_Pisa[spiders.selected.pisa[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Salt") {
				boot_Salt <- dataset_Salt[spiders.selected.salt[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Ther") {
				boot_Ther <- dataset_Ther[spiders.selected.ther[j,],]}

			if(Num_Per_Fam_Restricted$Family_Short %contains% "Thom") {
				boot_Thom <- dataset_Thom[spiders.selected.thom[j,],]}

			## Combine All Bootstrapped Datasets
			for (i in 1:nrow(Num_Per_Fam_Restricted)){
				parts <- get(paste0("boot_",Num_Per_Fam_Restricted$Family_Short[i]))
				boot <- rbind(boot,parts)
				}

			## Sort Dataset by Family_Short and also renames rows
			for (i in colnames(boot)){
				sql1 <- fn$identity("
					Select * From boot
					Order by Family_Short asc ")
				boot <- sqldf(c(sql1))
				}

			## Removes ID Columns
			boot_xid <- boot[2:24]

			## Cleanup
			if(exists("boot_Agel") == TRUE){rm(boot_Agel)}
			if(exists("boot_Cori") == TRUE){rm(boot_Cori)}
			if(exists("boot_Dict") == TRUE){rm(boot_Dict)}
			if(exists("boot_Gnap") == TRUE){rm(boot_Gnap)}
			if(exists("boot_Hahn") == TRUE){rm(boot_Hahn)}
			if(exists("boot_Liny") == TRUE){rm(boot_Liny)}
			if(exists("boot_Lyco") == TRUE){rm(boot_Lyco)}
			if(exists("boot_Pisa") == TRUE){rm(boot_Pisa)}
			if(exists("boot_Salt") == TRUE){rm(boot_Salt)}
			if(exists("boot_Ther") == TRUE){rm(boot_Ther)}
			if(exists("boot_Thom") == TRUE){rm(boot_Thom)}
			rm(parts)

		## Group Spiders Together to Generate Quantitative Data
			for (i in colnames(boot_xid)){
				sql1 <- fn$identity("
					Select Family_Short,
						sum(Agel) as Agel,
						sum(Cori) as Cori,
						sum(Dict) as Dict,
						sum(Gnap) as Gnap,
						sum(Hahn) as Hahn,
						sum(Liny) as Liny,
						sum(Lyco) as Lyco,
						sum(Pisa) as Pisa,
						sum(Salt) as Salt,
						sum(Ther) as Ther,
						sum(Thom) as Thom,
						sum(Arac) as Arac,
						sum(Cole) as Cole,
						sum(Coll) as Coll,
						sum(Derm) as Derm,
						sum(Dip) as Dip,
						sum(Gryl) as Gryl,
						sum(Hyme) as Hyme,
						sum(Isop) as Isop,
						sum(Lepi) as Lepi,
						sum(Opil) as Opil,
						sum(Pseu) as Pseu,
						count(Family_Short) as Total_Spiders
					From boot_xid
					Group by Family_Short ")
				boot_grouped <- sqldf(c(sql1))			## Grouped Dataset
					}
	
			Quantitative_frame <- boot_grouped[2:24]			## Quantitative Dataset as Data Frame
			row.names(Quantitative_frame) <- t(boot_grouped[1])	## Set the Row Names as the Family Names

		## Generate Indicator Dataset without X's where RowName = ColName (also removes total number of spiders column)
			Indicator_frame <- Quantitative_frame[1:22]
	
			## Use SQL to Change All Values > 0 to 1
				for (i in colnames(Indicator_frame)) {
					sql1 <- fn$identity("
						Update Indicator_frame
						Set $i = '$Link.yes'
						Where $i > '$Link.no' ")
					sql2 <- "Select * From main.Indicator_frame"
					Indicator_frame <- sqldf(c(sql1, sql2))		## Indicator Dataset as Data Frame
					}
	
				## Restore Row Names
				row.names(Indicator_frame) <- row.names(Quantitative_frame)	
	
		## Generate Bipartite Dataset where Predators = columns, Prey = Rows; Remove Total Spiders Column
			Quantitative_Bipartite <- t(Quantitative_frame[1:22])

		## Cleanup
			rm(boot)
			rm(boot_xid)

		## Calculate Directed Connectance From Indicator Data = TotalLinks / TotalPossibleLinks
			## TotalLinks <- sum(rowSums(as.matrix(Indicator_frame)))
			## TotalPossibleLinks <- nrow(as.matrix(Indicator_frame)) * ncol(as.matrix(Indicator_frame)) - nrow(as.matrix(Indicator_frame))	## Remove Cannibalism

			results.Connectance_directed[,j] <- sum(rowSums(as.matrix(Indicator_frame))) / ( nrow(as.matrix(Indicator_frame)) * ncol(as.matrix(Indicator_frame)) - nrow(as.matrix(Indicator_frame)) )

		## Interaction Strength From Quantitative Data
			## Generate Column Sums for Each Column
			## quant.colSum <- colSums(as.matrix(Quantitative_frame))
	
			## Calculate Prey Eaten as the Sum of the Number of Prey Item Columns
			## Prey_Eaten <- .rowSums(quant.colSum, m=1,n=22)
	
			## Calculate Total Number of Spiders for the Dataset by Taking Column Sum from the Last Column
			## Num_Spiders <- as.numeric(quant.colSum["Total_Spiders"])
	
			## Calculate Interaction Strength = Prey_Eaten / Num_Spiders
			results.Interaction_Strength[,j] <- .rowSums(colSums(as.matrix(Quantitative_frame)), m=1,n=22) / as.numeric(colSums(as.matrix(Quantitative_frame))["Total_Spiders"])

		## IntraGuildPredation Directed Connectance
			IGP_Conn_Pred_directed <- as.matrix(Indicator_frame)[,1:numberpredatoritems]
			IGP_Conn_Prey_directed <- as.matrix(Indicator_frame)[,(numberpredatoritems+1):(numberpredatoritems+numberpreyitems)]

			## Predators; Connectance = TotalLinks.pred / TotalPossibleLinks.pred
				## RowSum.Link.pred <- rowSums(IGP_Conn_Pred_directed)
				## TotalLinks.pred <- sum(RowSum.Link.pred)
				## TotalPossibleLinks.pred <- nrow(IGP_Conn_Pred_directed) * ncol(IGP_Conn_Pred_directed) - nrow(IGP_Conn_Pred_directed)	## Remove Cannibalism
			
				results.IGP_Conn_Pred_directed[,j] <- sum(rowSums(IGP_Conn_Pred_directed)) / ( nrow(IGP_Conn_Pred_directed) * ncol(IGP_Conn_Pred_directed) - nrow(IGP_Conn_Pred_directed) )

			## Prey; Connectance = TotalLinks.prey / TotalPossibleLinks.prey
				## RowSum.Link.prey <- rowSums(IGP_Conn_Prey_directed)
				## TotalLinks.prey <- sum(RowSum.Link.prey)
				## TotalPossibleLinks.prey <- nrow(IGP_Conn_Prey_directed) * ncol(IGP_Conn_Prey_directed)	## No Cannibalism present because dealing with strictly prey
			
				results.IGP_Conn_Prey_directed[,j] <- sum(rowSums(IGP_Conn_Prey_directed)) / ( nrow(IGP_Conn_Prey_directed) * ncol(IGP_Conn_Prey_directed) )

			## Cleanup
				rm(IGP_Conn_Pred_directed)
				rm(IGP_Conn_Prey_directed)

		## IntraGuildPredation Interaction Strength
			IGP_IS_Pred <- Quantitative_frame[c(1:numberpredatoritems,(numberpredatoritems+numberpreyitems+1))]
			IGP_IS_Prey <- Quantitative_frame[c((numberpredatoritems+1):(numberpredatoritems+numberpreyitems+1))]

			## Predators
				## Generate Column Sums for Each Column
				## quant.colSum.pred <- colSums(IGP_IS_Pred)
	
				## Calculate Prey Eaten as the Sum of the Number of Prey Item Columns
				## Prey_Eaten.pred <- .rowSums(quant.colSum.pred, m=1,n=numberpredatoritems)
		
				## Calculate Total Number of Spiders for the Dataset by Taking Column Sum from the Last Column
				## Num_Spiders.pred <- as.numeric(quant.colSum.pred[numberpredatoritems+1])
		
				## Calculate Interaction Strength = Prey_Eaten.pred / Num_Spiders.pred
				results.IGP_IS_Pred[,j] <- .rowSums(colSums(IGP_IS_Pred), m=1,n=numberpredatoritems) / as.numeric(colSums(IGP_IS_Pred)[numberpredatoritems+1])
	
			## Prey
				## Generate Column Sums for Each Column
				## quant.colSum.prey <- colSums(IGP_IS_Prey)
	
				## Calculate Prey Eaten as the Sum of the Number of Prey Item Columns
				## Prey_Eaten.prey <- .rowSums(quant.colSum.prey, m=1,n=numberpreyitems)
		
				## Calculate Total Number of Spiders for the Dataset by Taking Column Sum from the Last Column
				## Num_Spiders.prey <- as.numeric(quant.colSum.prey[numberpreyitems+1])
		
				## Calculate Interaction Strength = Prey_Eaten.prey / Num_Spiders.prey
				results.IGP_IS_Prey[,j] <- .rowSums(colSums(IGP_IS_Prey), m=1,n=numberpreyitems) / as.numeric(colSums(IGP_IS_Prey)[numberpreyitems+1])

			## Cleanup
				rm(IGP_IS_Pred)
				rm(IGP_IS_Prey)

		## Compartmentalization from Bipartite
			results.Compart_Bipartite[,j] <- compart(Quantitative_Bipartite)$n.compart

		## Interaction Evenness Bipartite
			## logbase e; intereven = sum
			results.Interaction_Evenness_e_sum[,j] <- as.numeric(
				networklevel(Quantitative_Bipartite,
					index="interaction evenness", 
					level="both",
					weighted=TRUE,
					ISAmethod="Bluethgen",
					SAmethod = "Bluethgen",
					extinctmethod = "random",
					nrep = 100,
					CCfun=median,
					dist="horn",
					normalise=TRUE,
					empty.web=TRUE,
					logbase="e",
					intereven="sum",
					H2_integer=TRUE,
					fcweighted=TRUE,
					fcdist="euclidean",
					legacy=FALSE))

			## logbase 2; intereven = sum
			results.Interaction_Evenness_2_sum[,j] <- as.numeric(
				networklevel(Quantitative_Bipartite,
					index="interaction evenness", 
					level="both",
					weighted=TRUE,
					ISAmethod="Bluethgen",
					SAmethod = "Bluethgen",
					extinctmethod = "random",
					nrep = 100,
					CCfun=median,
					dist="horn",
					normalise=TRUE,
					empty.web=TRUE,
					logbase="2",
					intereven="sum",
					H2_integer=TRUE,
					fcweighted=TRUE,
					fcdist="euclidean",
					legacy=FALSE))

			## logbase e; intereven = prod
			results.Interaction_Evenness_e_prod[,j] <- as.numeric(
				networklevel(Quantitative_Bipartite,
					index="interaction evenness", 
					level="both",
					weighted=TRUE,
					ISAmethod="Bluethgen",
					SAmethod = "Bluethgen",
					extinctmethod = "random",
					nrep = 100,
					CCfun=median,
					dist="horn",
					normalise=TRUE,
					empty.web=TRUE,
					logbase="e",
					intereven="prod",
					H2_integer=TRUE,
					fcweighted=TRUE,
					fcdist="euclidean",
					legacy=FALSE))

			## logbase 2; intereven = prod
			results.Interaction_Evenness_2_prod[,j] <- as.numeric(
				networklevel(Quantitative_Bipartite,
					index="interaction evenness", 
					level="both",
					weighted=TRUE,
					ISAmethod="Bluethgen",
					SAmethod = "Bluethgen",
					extinctmethod = "random",
					nrep = 100,
					CCfun=median,
					dist="horn",
					normalise=TRUE,
					empty.web=TRUE,
					logbase="2",
					intereven="prod",
					H2_integer=TRUE,
					fcweighted=TRUE,
					fcdist="euclidean",
					legacy=FALSE))

		## Specialization H2 and H2 Uncorrelated Prey (Pred = col, Prey = row)
			H2fun_results <- H2fun(web = Quantitative_Bipartite, H2_integer = TRUE)
			results.Specialization_h2[,j] <- as.numeric(H2fun_results[1])
			results.Specialization_h2_uncorr[,j] <- as.numeric(H2fun_results[4])

			## Cleanup
			rm(H2fun_results)

		## Specialization H2 and H2 Uncorrelated Predators (Prey = col, Pred = row)
			H2fun_results_pred <- H2fun(web = t(Quantitative_Bipartite), H2_integer = TRUE)
			results.Specialization_pred_h2[,j] <- as.numeric(H2fun_results_pred[1])
			results.Specialization_pred_h2_uncorr[,j] <- as.numeric(H2fun_results_pred[4])

			## Cleanup
			rm(H2fun_results_pred)

		## Specialization dprime Prey (Pred = col, Prey = row)
			abundances <- t(as.matrix(Quantitative_frame)[,(numberpredatoritems+numberpreyitems+1):(numberpredatoritems+numberpreyitems+1)])
			## dfun_results <- dfun(web = Quantitative_Bipartite, abuns = abundances)
			results.Specialization_dprime[,j] <- t(t(dfun(web = Quantitative_Bipartite, abuns = abundances)$dprime))

			## Cleanup
			rm(abundances)

		## Specialization dprime Pred (dprime transpose; Prey = col, Pred = row)
			if(exists("results.Specialization_pred_dprime") == FALSE){
				## Specialization predators dprime Results Matrix Creation
				results.Specialization_pred_dprime <- matrix(, nrow = nrow(as.matrix(Quantitative_frame)), ncol = B, byrow = TRUE)
				rownames(results.Specialization_pred_dprime) <- rownames(as.matrix(Quantitative_frame))
				}

			dfun_results_pred <- dfun(web = t(Quantitative_Bipartite), abuns = NULL)
			if(nrow(t(t(dfun_results_pred$dprime))) == nrow(results.Specialization_pred_dprime)){
				results.Specialization_pred_dprime[,j] <- t(t(dfun_results_pred$dprime))
				
				## Cleanup
				rm(dfun_results_pred)
				}else{

				## When the results don't have all of the predators, add them in and sort them
				temp5 <- t(t(dfun_results_pred$dprime))

				if((rownames(temp5) %contains% "Agel" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Agel" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Agel")
					}

				if((rownames(temp5) %contains% "Cori" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Cori" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Cori")
					}

				if((rownames(temp5) %contains% "Dict" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Dict" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Dict")
					}

				if((rownames(temp5) %contains% "Gnap" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Gnap" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Gnap")
					}

				if((rownames(temp5) %contains% "Hahn" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Hahn" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Hahn")
					}

				if((rownames(temp5) %contains% "Liny" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Liny" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Liny")
					}

				if((rownames(temp5) %contains% "Lyco" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Lyco" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Lyco")
					}

				if((rownames(temp5) %contains% "Pisa" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Pisa" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Pisa")
					}

				if((rownames(temp5) %contains% "Salt" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Salt" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Salt")
					}

				if((rownames(temp5) %contains% "Ther" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Ther" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Ther")
					}

				if((rownames(temp5) %contains% "Thom" == FALSE) & (Num_Per_Fam_Restricted$Family_Short %contains% "Thom" == TRUE)){
					temp5_row_names_temp <- row.names(temp5)
					temp5 <- rbind(temp5, NaN)
					row.names(temp5) <- c(temp5_row_names_temp, "Thom")
					}

				## Sort temp5
				temp6 <- t(t(temp5[order(rownames(temp5)),]))

				## Store results
				results.Specialization_pred_dprime[,j] <- temp6

				## Cleanup
				rm(dfun_results_pred)
				rm(temp5)
				rm(temp5_row_names_temp)
				rm(temp6)
				}

		## Undirected Indicator Dataset with All Predators (For use with Connectance and Compartmentalization Pimm & Lawton)
			temp1 <- Indicator_frame
	
			## Add Predators that Don't currently exist in the Directed Dataset
			if(rownames(temp1) %contains% "Agel" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Agel")
				}
	
			if(rownames(temp1) %contains% "Cori" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Cori")
				}
	
			if(rownames(temp1) %contains% "Dict" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Dict")
				}
	
			if(rownames(temp1) %contains% "Gnap" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Gnap")
				}
	
			if(rownames(temp1) %contains% "Hahn" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Hahn")
				}
	
			if(rownames(temp1) %contains% "Liny" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Liny")
				}
	
			if(rownames(temp1) %contains% "Lyco" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Lyco")
				}
	
			if(rownames(temp1) %contains% "Pisa" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Pisa")
				}
	
			if(rownames(temp1) %contains% "Salt" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Salt")
				}
	
			if(rownames(temp1) %contains% "Ther" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Ther")
				}
	
			if(rownames(temp1) %contains% "Thom" == FALSE){
				temp1_row_names_temp <- row.names(temp1)
				temp1 <- rbind(temp1,vec0)
				row.names(temp1) <- c(temp1_row_names_temp,"Thom")
				}
	
			## Sort Dataset by Family and also renames rows
			if(nrow(Num_Per_Fam_Restricted) < 11){
			temp1_row_names <- row.names(temp1)			## Generates Final List of Row Names
			temp1["Family_Short"] <- temp1_row_names		## Stores List of Row Names in Family_Short Column
	
			for (i in colnames(temp1)){
				sql1 <- fn$identity("
					Select * From temp1
					Order by Family_Short asc ")
				temp1_sorted <- sqldf(c(sql1))		## Sorts temp1
				}
	
			row.names(temp1_sorted) <- t(temp1_sorted["Family_Short"])	## Gives Row Names to Sorted Temp1
	
			temp2 <- temp1_sorted[1:22]				## Gets rid of Family_Short Column in Temp1
	
	
			## Verify that Row Names Match First Predator Column Names
			if(is.logical(rownames(temp2) == colnames(temp2[1:numberpredatoritems])) == FALSE){break}
			}else{temp2 <- Indicator_frame}

			## Create Undirected Indicator Matrix (for use with Pimm & Lawton, and Undirected Connectance)
			Indicator_undirected_matrix <- matrix(data = 0, nrow = ncol(temp2), ncol = ncol(temp2))
	
			## Update Undirected Indicator Matrix so that all links are shown
			## Predators = rows, Prey = columns
			## Predators = k, Prey = l
			for(k in 1:nrow(temp2)){
				for(l in 1:ncol(temp2)){
					if(temp2[k,l] == 1){(Indicator_undirected_matrix[k,l] <- 1) & (Indicator_undirected_matrix[l,k] <- 1)}
				}
			}
	
			## Test to verify if symmetric or not.  If not symmetric, then process stops because of error.
				if(isSymmetric(Indicator_undirected_matrix) == FALSE){break}
	
			## Add Back Row Names
			rownames(Indicator_undirected_matrix) <- colnames(temp2)
			colnames(Indicator_undirected_matrix) <- colnames(temp2)

			## Cleanup
			if(nrow(Num_Per_Fam_Restricted) < 11){
			rm(temp1_sorted)
			rm(temp1_row_names_temp)
			rm(temp1_row_names)
			}else{rm(temp1)}

		## Compartmentalization Pimm & Lawton - All Spiders (Transpose not used because matrix is symmetric)
			## Make a Row Link Matrix and Break into Individual Rows
			pimm.all.link.row.matrix <- matrix(data = NA, nrow = nrow(Indicator_undirected_matrix), ncol = ncol(Indicator_undirected_matrix))
			for(i in 1:nrow(Indicator_undirected_matrix)){
				for( k in 1:ncol(Indicator_undirected_matrix)){
					if(Indicator_undirected_matrix[i,k] == 1) pimm.all.link.row.matrix [i,k] <- k
		 			}
				}
	
			## Make a Column Link Matrix and Break into Individual Columns
			pimm.all.link.col.matrix <- matrix(data = NA, nrow = nrow(Indicator_undirected_matrix), ncol = ncol(Indicator_undirected_matrix))
			for(i in 1:nrow(Indicator_undirected_matrix)){
				for(k in 1:ncol(Indicator_undirected_matrix)){
					if(Indicator_undirected_matrix[i,k] == 1) pimm.all.link.col.matrix[i,k] <- i
					}
				}

			## Calculating Numerator
			## Number of species with which both i = col, j = row interact
			pimm.all.numerator <- matrix(data = NA, nrow = nrow(Indicator_undirected_matrix), ncol = ncol(Indicator_undirected_matrix))
			for(k in 1:nrow(pimm.all.link.row.matrix)){
				for(i in 1:ncol(pimm.all.link.col.matrix)){
					pimm.all.numerator[k,i] <- length(intersect(pimm.all.link.row.matrix[k,], pimm.all.link.col.matrix[,i])) - 1 ## Subtract 1 because it counts all the NA's as a value, which we want to omit
					}
				}
		
			## Calculating Denominator
			## Number of species which which either i = col, j = row interact
			## Take the number of links minus the duplicates, which would be found in the numerator matrix
			pimm.all.denominator <- matrix(data = NA, nrow = nrow(Indicator_undirected_matrix), ncol = ncol(Indicator_undirected_matrix))
			for(k in 1:nrow(pimm.all.link.row.matrix)){
				for(i in 1:ncol(pimm.all.link.col.matrix)){
					pimm.all.denominator[k,i] <- length(union(pimm.all.link.row.matrix[k,], pimm.all.link.col.matrix[,i])) - 1 ## Subtract 1 because it counts all the NA's as a value, which we want to omit
					}
				}
		
			## Calculating c[j,i]
			pimm.all.stat.c <- pimm.all.numerator/pimm.all.denominator
			for(i in 1:nrow(pimm.all.stat.c)){
				pimm.all.stat.c[i,i] <- 0
				}
	
			## Updating c[j,i] to 0 if there are NAs
			for(i in 1:nrow(pimm.all.stat.c)){
				for(k in 1:ncol(pimm.all.stat.c)){
					if(pimm.all.stat.c[i,k] == "NaN"){ pimm.all.stat.c[i,k] <- 0 }
				}
			}
	
			## Calculating Compartmentalization Statistic
			results.Compartmentalization_pimm_all[,j] <- sum(pimm.all.stat.c) / ( nrow(pimm.all.stat.c) * (nrow(pimm.all.stat.c) - 1) )
	
			## Cleanup
			rm(pimm.all.link.row.matrix)
			rm(pimm.all.link.col.matrix)
			rm(pimm.all.numerator)
			rm(pimm.all.denominator)
			rm(pimm.all.stat.c)


		## Compartmentalization Pimm & Lawton - Predators Only (Transpose not used because matrix is symmetric)
			temp2.short.pred <- temp2[,1:numberpredatoritems]
	
			## Create Compartmentalization Matrix
			pimm.matrix.pred <- matrix(data = 0, nrow = nrow(temp2.short.pred), ncol = numberpredatoritems)
	
			## Update Compartmentalization Matrix so that all links are shown
			## Predators = rows, Prey = columns
			## Predators = k, Prey = l
			for(k in 1:nrow(temp2.short.pred)){
				for(l in 1:ncol(temp2.short.pred)){
					if(temp2.short.pred[k,l] == 1){(pimm.matrix.pred[k,l] <- 1) & (pimm.matrix.pred[l,k] <- 1)}
				}
			}
	
			## Test to verify if symmetric or not.  If not symmetric, then process stops because of error.
			if(isSymmetric(pimm.matrix.pred) == FALSE){break}
	
			## Add Back Row Names
			rownames(pimm.matrix.pred) <- rownames(temp2.short.pred)
			colnames(pimm.matrix.pred) <- colnames(temp2.short.pred)
		
			## Make a Row Link Matrix and Break into Individual Rows
			pimm.link.row.matrix.pred <- matrix(data = NA, nrow = nrow(pimm.matrix.pred), ncol = ncol(pimm.matrix.pred))
			for(i in 1:nrow(pimm.matrix.pred)){
				for(k in 1:ncol(pimm.matrix.pred)){
					if(pimm.matrix.pred[i,k] == 1) pimm.link.row.matrix.pred[i,k] <- k
		 			}
				}
	
			## Make a Column Link Matrix and Break into Individual Columns
			pimm.link.col.matrix.pred <- matrix(data = NA, nrow = nrow(pimm.matrix.pred), ncol = ncol(pimm.matrix.pred))
			for(i in 1:nrow(pimm.matrix.pred)){
				for(k in 1:ncol(pimm.matrix.pred)){
					if(pimm.matrix.pred[i,k] == 1) pimm.link.col.matrix.pred[i,k] <- i
					}
				}
	
			## Calculating Numerator
			## Number of species with which both i = col, j = row interact
			pimm.numerator.pred <- matrix(data = NA, nrow = nrow(pimm.matrix.pred), ncol = ncol(pimm.matrix.pred))
			for(k in 1:nrow(pimm.link.row.matrix.pred)){
				for(i in 1:ncol(pimm.link.col.matrix.pred)){
					pimm.numerator.pred[k,i] <- length(intersect(pimm.link.row.matrix.pred[k,], pimm.link.col.matrix.pred[,i])) - 1 ## Subtract 1 because it counts all the NA's as a value, which we want to omit
					}
				}
	
			## Calculating Denominator
			## Number of species which which either i = col, j = row interact
			## Take the number of links minus the duplicates, which would be found in the numerator matrix
			pimm.denominator.pred <- matrix(data = NA, nrow = nrow(pimm.matrix.pred), ncol = ncol(pimm.matrix.pred))
			for(k in 1:nrow(pimm.link.row.matrix.pred)){
				for(i in 1:ncol(pimm.link.col.matrix.pred)){
					pimm.denominator.pred[k,i] <- length(union(pimm.link.row.matrix.pred[k,], pimm.link.col.matrix.pred[,i])) - 1 ## Subtract 1 because it counts all the NA's as a value, which we want to omit
					}
				}
	
			## Calculating c[j,i]
			pimm.stat.c.pred <- pimm.numerator.pred/pimm.denominator.pred
			for(i in 1:nrow(pimm.stat.c.pred)){
				pimm.stat.c.pred[i,i] <- 0
				}
	
			## Updating c[j,i] to 0 if there are NAs
			for(i in 1:nrow(pimm.stat.c.pred)){
				for(k in 1:ncol(pimm.stat.c.pred)){
					if(pimm.stat.c.pred[i,k] == "NaN"){ pimm.stat.c.pred[i,k] <- 0 }
				}
			}
	
			## Calculating Compartmentalization Statistic
			results.Compartmentalization_pimm_pred[,j] <- sum(pimm.stat.c.pred) / ( nrow(pimm.stat.c.pred) * (nrow(pimm.stat.c.pred) - 1) )

			## Cleanup
			rm(temp2.short.pred)
			rm(temp2)
			rm(pimm.link.row.matrix.pred)
			rm(pimm.link.col.matrix.pred)
			rm(pimm.numerator.pred)
			rm(pimm.denominator.pred)
			rm(pimm.stat.c.pred)

		## Calculate Undirected Connectance From Indicator Data = L_all / ( S_all * (S_all - 1) )
			## L_all <- sum(colSums(Indicator_undirected_matrix))
			## S_all <- numberpredatoritems + numberpreyitems
			results.Connectance_undirected[,j] <- sum(colSums(Indicator_undirected_matrix)) / ( (numberpredatoritems + numberpreyitems) * ((numberpredatoritems + numberpreyitems) - 1) )

		## Calculated Undirected Connectance (including cannibalism) From Indicator Data = L_all / ( S_all^2 )
			## L_all <- sum(colSums(Indicator_undirected_matrix))
			##S_all <- numberpredatoritems + numberpreyitems
			results.Connectance_wcannibalism[,j] <- sum(colSums(Indicator_undirected_matrix)) / ( (numberpredatoritems + numberpreyitems)^2 )

		## IntraGuildPredation Undirected Connectance Predators Only (Prey doesn't make sense to do here) = L_IGP_pred / ( S_IGP_pred * (S_IGP_pred - 1) )
			IGP_Conn_Pred_undirected <- Indicator_undirected_matrix[1:numberpredatoritems,1:numberpredatoritems]
	
			## L_IGP_pred <- sum(colSums(IGP_Conn_Pred_undirected))
			## S_IGP_pred <- numberpredatoritems
			results.IGP_Conn_Pred_undirected[,j] <- sum(colSums(IGP_Conn_Pred_undirected)) / ( numberpredatoritems * (numberpredatoritems - 1) )
	
			## Cleanup
			rm(IGP_Conn_Pred_undirected)
	
		## Contingency Table
		if(data_type == "ALL"){
			## Splitting Datasets
				for(i in colnames(boot_grouped)){
					sql1 <- fn$identity("
						Select Family_Short,
							sum(Agel) as Agel,
							sum(Cori) as Cori,
							sum(Dict) as Dict,
							sum(Gnap) as Gnap,
							sum(Hahn) as Hahn,
							sum(Liny) as Liny,
							sum(Lyco) as Lyco,
							sum(Pisa) as Pisa,
							sum(Salt) as Salt,
							sum(Ther) as Ther,
							sum(Thom) as Thom
						From boot_grouped
						Where Family_Short in ('Cori', 'Gnap', 'Lyco', 'Pisa', 'Salt', 'Thom')
						Group by Family_Short ")
					sql2 <- fn$identity("
						Select Family_Short,
							sum(Arac) as Arac,
							sum(Cole) as Cole,
							sum(Coll) as Coll,
							sum(Derm) as Derm,
							sum(Dip) as Dip,
							sum(Gryl) as Gryl,
							sum(Hyme) as Hyme,
							sum(Isop) as Isop,
							sum(Lepi) as Lepi,
							sum(Opil) as Opil,
							sum(Pseu) as Pseu
						From boot_grouped
						Where Family_Short in ('Cori', 'Gnap', 'Lyco', 'Pisa', 'Salt', 'Thom')
						Group by Family_Short ")
					sql3 <- fn$identity("
						Select Family_Short,
							sum(Agel) as Agel,
							sum(Cori) as Cori,
							sum(Dict) as Dict,
							sum(Gnap) as Gnap,
							sum(Hahn) as Hahn,
							sum(Liny) as Liny,
							sum(Lyco) as Lyco,
							sum(Pisa) as Pisa,
							sum(Salt) as Salt,
							sum(Ther) as Ther,
							sum(Thom) as Thom
						From boot_grouped
						Where Family_Short in ('Agel', 'Dict', 'Hahn', 'Liny', 'Ther')
						Group by Family_Short ")
					sql4 <- fn$identity("
						Select Family_Short,
							sum(Arac) as Arac,
							sum(Cole) as Cole,
							sum(Coll) as Coll,
							sum(Derm) as Derm,
							sum(Dip) as Dip,
							sum(Gryl) as Gryl,
							sum(Hyme) as Hyme,
							sum(Isop) as Isop,
							sum(Lepi) as Lepi,
							sum(Opil) as Opil,
							sum(Pseu) as Pseu
						From boot_grouped
						Where Family_Short in ('Agel', 'Dict', 'Hahn', 'Liny', 'Ther')
						Group by Family_Short ")
					dataset_curs_spiders <- sqldf(c(sql1))
					dataset_curs_nonspiders <- sqldf(c(sql2))
					dataset_web_spiders <- sqldf(c(sql3))
					dataset_web_nonspiders <- sqldf(c(sql4))
					}
	
					dataset_curs_spiders2 <- dataset_curs_spiders[,2:(numberpredatoritems+1)]
					dataset_curs_nonspiders2 <- dataset_curs_nonspiders[,2:(numberpreyitems+1)]
					dataset_web_spiders2 <- dataset_web_spiders[,2:(numberpredatoritems+1)]
					dataset_web_nonspiders2 <- dataset_web_nonspiders[,2:(numberpreyitems+1)]
	
					row.names(dataset_curs_spiders2) <- t(dataset_curs_spiders["Family_Short"])
					row.names(dataset_curs_nonspiders2) <- t(dataset_curs_nonspiders["Family_Short"])
					row.names(dataset_web_spiders2) <- t(dataset_web_spiders["Family_Short"])
					row.names(dataset_web_nonspiders2) <- t(dataset_web_nonspiders["Family_Short"])
	
	
			## Changing Split Datasets to Indicators
				indicator_curs_spiders <- dataset_curs_spiders2
				indicator_curs_nonspiders <- dataset_curs_nonspiders2
				indicator_web_spiders <- dataset_web_spiders2
				indicator_web_nonspiders <- dataset_web_nonspiders2
	
				indicator_curs_spiders_rownames <- row.names(indicator_curs_spiders)
				indicator_curs_nonspiders_rownames <- row.names(indicator_curs_nonspiders)
				indicator_web_spiders_rownames <- row.names(indicator_web_spiders)
				indicator_web_nonspiders_rownames <- row.names(indicator_web_nonspiders)
	
	
				for (i in colnames(indicator_curs_spiders)) {
					sql1 <- fn$identity("
						Update indicator_curs_spiders
						Set $i = '$Link.yes'
						Where $i > '$Link.no' ")
					sql2 <- "Select * From indicator_curs_spiders"
					indicator_curs_spiders <- sqldf(c(sql1,sql2))
					}
	
				for (i in colnames(indicator_curs_nonspiders)) {
					sql1 <- fn$identity("
						Update indicator_curs_nonspiders
						Set $i = '$Link.yes'
						Where $i > '$Link.no' ")
					sql2 <- "Select * From indicator_curs_nonspiders"
					indicator_curs_nonspiders <- sqldf(c(sql1,sql2))
					}
	
				for (i in colnames(indicator_web_spiders)) {
					sql1 <- fn$identity("
						Update indicator_web_spiders
						Set $i = '$Link.yes'
						Where $i > '$Link.no' ")
					sql2 <- "Select * From indicator_web_spiders"
					indicator_web_spiders <- sqldf(c(sql1,sql2))
					}
	
				for (i in colnames(indicator_web_nonspiders)) {
					sql1 <- fn$identity("
						Update indicator_web_nonspiders
						Set $i = '$Link.yes'
						Where $i > '$Link.no' ")
					sql2 <- "Select * From indicator_web_nonspiders"
					indicator_web_nonspiders <- sqldf(c(sql1,sql2))
					}
	
				row.names(indicator_curs_spiders) <- indicator_curs_spiders_rownames
				row.names(indicator_curs_nonspiders) <- indicator_curs_nonspiders_rownames
				row.names(indicator_web_spiders) <- indicator_web_spiders_rownames
				row.names(indicator_web_nonspiders) <- indicator_web_nonspiders_rownames
	
	
			## Count of Spiders Quant
				quant_count_curs_spiders <- sum(dataset_curs_spiders2)
				quant_count_curs_nonspiders <- sum(dataset_curs_nonspiders2)
				quant_count_web_spiders <- sum(dataset_web_spiders2)
				quant_count_web_nonspiders <- sum(dataset_web_nonspiders2)
	
	
			## Count of Spiders Indicator
				ind_count_curs_spiders <- sum(indicator_curs_spiders)
				ind_count_curs_nonspiders <- sum(indicator_curs_nonspiders)
				ind_count_web_spiders <- sum(indicator_web_spiders)
				ind_count_web_nonspiders <- sum(indicator_web_nonspiders)
	
			## Cleanup
				rm(dataset_curs_spiders)
				rm(dataset_curs_nonspiders)
				rm(dataset_web_spiders)
				rm(dataset_web_nonspiders)
				rm(dataset_curs_spiders2)
				rm(dataset_curs_nonspiders2)
				rm(dataset_web_spiders2)
				rm(dataset_web_nonspiders2)
				rm(indicator_curs_spiders)
				rm(indicator_curs_nonspiders)
				rm(indicator_web_spiders)
				rm(indicator_web_nonspiders)
				rm(indicator_curs_spiders_rownames)
				rm(indicator_curs_nonspiders_rownames)
				rm(indicator_web_spiders_rownames)
				rm(indicator_web_nonspiders_rownames)
	
	
			## Contingency Table Quant
				quant_cont_table <- matrix(c(
						quant_count_curs_spiders, quant_count_curs_nonspiders,
						quant_count_web_spiders, quant_count_web_nonspiders),
					nrow = 2, ncol = 2,
					byrow = TRUE,
					dimnames = list(c("Cursorial","Web Building"),
								c("Spiders", "Non-Spiders")))
	
				quant_fisher_central <- fisher.exact(quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "central")
				quant_fisher_minlike <- fisher.exact(quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "minlike")
				quant_blaker <- blaker.exact(quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05)
	
				results.quant_fisher_central_oddsratio[,j] <- as.numeric(quant_fisher_central$estimate)
				results.quant_fisher_central_pval[,j] <- quant_fisher_central$p.value
				results.quant_fisher_central_confint[j,1] <- quant_fisher_central$conf.int[1]
				results.quant_fisher_central_confint[j,2] <- quant_fisher_central$conf.int[2]
	
				results.quant_fisher_minlike_oddsratio[,j] <- as.numeric(quant_fisher_minlike$estimate)
				results.quant_fisher_minlike_pval[,j] <- quant_fisher_minlike$p.value
				results.quant_fisher_minlike_confint[j,1] <- quant_fisher_minlike$conf.int[1]
				results.quant_fisher_minlike_confint[j,2] <- quant_fisher_minlike$conf.int[2]
	
				results.quant_blaker_oddsratio[,j] <- as.numeric(quant_blaker$estimate)
				results.quant_blaker_pval[,j] <- quant_blaker$p.value
				results.quant_blaker_confint[j,1] <- quant_blaker$conf.int[1]
				results.quant_blaker_confint[j,2] <- quant_blaker$conf.int[2]
	
				## Cleanup
				rm(quant_fisher_central)
				rm(quant_fisher_minlike)
				rm(quant_blaker)
				rm(quant_count_curs_spiders)
				rm(quant_count_curs_nonspiders)
				rm(quant_count_web_spiders)
				rm(quant_count_web_nonspiders)
				rm(quant_cont_table)
	
			## Contingency Table Indicator
				ind_cont_table <- matrix(c(
						ind_count_curs_spiders, ind_count_curs_nonspiders,
						ind_count_web_spiders, ind_count_web_nonspiders),
					nrow = 2, ncol = 2,
					byrow = TRUE,
					dimnames = list(c("Cursorial", "Web Building"),
								c("Spiders", "Non-Spiders")))
	
				ind_fisher_central <- fisher.exact(ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "central")
				ind_fisher_minlike <- fisher.exact(ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "minlike")
				ind_blaker <- blaker.exact(ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05)
	
				results.ind_fisher_central_oddsratio[,j] <- as.numeric(ind_fisher_central$estimate)
				results.ind_fisher_central_pval[,j] <- ind_fisher_central$p.value
				results.ind_fisher_central_confint[j,1] <- ind_fisher_central$conf.int[1]
				results.ind_fisher_central_confint[j,2] <- ind_fisher_central$conf.int[2]
	
				results.ind_fisher_minlike_oddsratio[,j] <- as.numeric(ind_fisher_minlike$estimate)
				results.ind_fisher_minlike_pval[,j] <- ind_fisher_minlike$p.value
				results.ind_fisher_minlike_confint[j,1] <- ind_fisher_minlike$conf.int[1]
				results.ind_fisher_minlike_confint[j,2] <- ind_fisher_minlike$conf.int[2]
	
				results.ind_blaker_oddsratio[,j] <- as.numeric(ind_blaker$estimate)
				results.ind_blaker_pval[,j] <- ind_blaker$p.value
				results.ind_blaker_confint[j,1] <- ind_blaker$conf.int[1]
				results.ind_blaker_confint[j,2] <- ind_blaker$conf.int[2]
	
				## Cleanup
				rm(ind_fisher_central)
				rm(ind_fisher_minlike)
				rm(ind_blaker)
				rm(ind_count_curs_spiders)
				rm(ind_count_curs_nonspiders)
				rm(ind_count_web_spiders)
				rm(ind_count_web_nonspiders)
				rm(ind_cont_table)
	
			} ## END OF CONTINGENCY TABLE
	
		## Progress Bar
			   Sys.sleep(0.1)
			   setTkProgressBar(pb_boot, j, title=paste( round(j/B*100, 0),"% done"))

		## Get Rid of All Excess Stored Results except for the last time
		if((j == B) == FALSE){	
			rm(boot_grouped)
			rm(Indicator_frame)
			rm(Indicator_undirected_matrix)
			rm(pimm.matrix.pred)
			rm(Quantitative_frame)
			rm(Quantitative_Bipartite)
			}

		} ## END OF BOOTSTRAP
close(pb_boot)
################
################
################
################
################
################ BOOTSTRAP RESULTS
## Sample Variances
#	var.Connectance <- var(results.Connectance)
#	var.Interaction_Strength <- var(results.Interaction_Strength)
#	var.IGP_Conn_Pred <- var(results.IGP_Conn_Pred)
#	var.IGP_Conn_Prey <- var(results.IGP_Conn_Prey)
#	var.IGP_IS_Pred <- var(results.IGP_IS_Pred)
#	var.IGP_IS_Prey <- var(results.IGP_IS_Prey)
#	var.Compart_Bipartite <- var(results.Compart_Bipartite)
#	var.Interaction_Evenness_e_sum <- var(results.Interaction_Evenness_e_sum)
#	var.Interaction_Evenness_2_sum <- var(results.Interaction_Evenness_2_sum)
#	var.Interaction_Evenness_e_prod <- var(results.Interaction_Evenness_e_prod)
#	var.Interaction_Evenness_2_prod <- var(results.Interaction_Evenness_2_prod)


## Sample Medians
	med.Connectance_wcannibalism <- median(results.Connectance_wcannibalism)
	med.Connectance_undirected <- median(results.Connectance_undirected)
	med.Connectance_directed <- median(results.Connectance_directed)
	med.Interaction_Strength <- median(results.Interaction_Strength)
	med.IGP_Conn_Pred_undirected <- median(results.IGP_Conn_Pred_undirected)
	med.IGP_Conn_Pred_directed <- median(results.IGP_Conn_Pred_directed)
	med.IGP_Conn_Prey_directed <- median(results.IGP_Conn_Prey_directed)
	med.IGP_IS_Pred <- median(results.IGP_IS_Pred)
	med.IGP_IS_Prey <- median(results.IGP_IS_Prey)
	med.Interaction_Evenness_e_sum <- median(results.Interaction_Evenness_e_sum)
	med.Interaction_Evenness_2_sum <- median(results.Interaction_Evenness_2_sum)
	med.Interaction_Evenness_e_prod <- median(results.Interaction_Evenness_e_prod)
	med.Interaction_Evenness_2_prod <- median(results.Interaction_Evenness_2_prod)
	med.Specialization_h2 <- median(results.Specialization_h2)
	med.Specialization_h2_uncorr <- median(results.Specialization_h2_uncorr)
	med.Specialization_pred_h2 <- median(results.Specialization_pred_h2)
	med.Specialization_pred_h2_uncorr <- median(results.Specialization_pred_h2_uncorr)
	med.Compart_Bipartite <- median(results.Compart_Bipartite)
	med.Compartmentalization_pimm_all <- median(results.Compartmentalization_pimm_all)
	med.Compartmentalization_pimm_pred <- median(results.Compartmentalization_pimm_pred)


## Percentile Bootstrap Confidence Interval
	CIp.Connectance_wcannibalism <- as.numeric(quantile(results.Connectance_wcannibalism, c(alpha / 2, 1 - alpha / 2)))
	CIp.Connectance_undirected <- as.numeric(quantile(results.Connectance_undirected, c(alpha / 2, 1 - alpha / 2)))
	CIp.Connectance_directed <- as.numeric(quantile(results.Connectance_directed, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Interaction_Strength <- as.numeric(quantile(results.Interaction_Strength, c(alpha / 2 , 1 - alpha / 2)))
	CIp.IGP_Conn_Pred_undirected <- as.numeric(quantile(results.IGP_Conn_Pred_undirected, c(alpha / 2 , 1 - alpha / 2)))
	CIp.IGP_Conn_Pred_directed <- as.numeric(quantile(results.IGP_Conn_Pred_directed, c(alpha / 2 , 1 - alpha / 2)))
	CIp.IGP_Conn_Prey_directed <- as.numeric(quantile(results.IGP_Conn_Prey_directed, c(alpha / 2 , 1 - alpha / 2)))
	CIp.IGP_IS_Pred <- as.numeric(quantile(results.IGP_IS_Pred, c(alpha / 2 , 1 - alpha / 2)))
	CIp.IGP_IS_Prey <- as.numeric(quantile(results.IGP_IS_Prey, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Interaction_Evenness_e_sum <- as.numeric(quantile(results.Interaction_Evenness_e_sum, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Interaction_Evenness_2_sum <- as.numeric(quantile(results.Interaction_Evenness_2_sum, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Interaction_Evenness_e_prod <- as.numeric(quantile(results.Interaction_Evenness_e_prod, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Interaction_Evenness_2_prod <- as.numeric(quantile(results.Interaction_Evenness_2_prod, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Specialization_h2 <- as.numeric(quantile(results.Specialization_h2, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Specialization_h2_uncorr <- as.numeric(quantile(results.Specialization_h2_uncorr, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Specialization_pred_h2 <- as.numeric(quantile(results.Specialization_pred_h2, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Specialization_pred_h2_uncorr <- as.numeric(quantile(results.Specialization_pred_h2_uncorr, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Compart_Bipartite <- as.numeric(quantile(results.Compart_Bipartite, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Compartmentalization_pimm_all <- as.numeric(quantile(results.Compartmentalization_pimm_all, c(alpha / 2 , 1 - alpha / 2)))
	CIp.Compartmentalization_pimm_pred <- as.numeric(quantile(results.Compartmentalization_pimm_pred, c(alpha / 2 , 1 - alpha / 2)))


## Bootstrap "T" Confidence Interval
	## This requires the standard deviation of the values from the original model, which we don't have.


## Confidence Intervals per row for Specialization dprime and dprime transpose
	CIp.Specialization_dprime <- rowQuantiles(results.Specialization_dprime, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, drop = FALSE)			## Na's are excluded first, do not drop singleton dimensions in the result
		rownames(CIp.Specialization_dprime) <- rownames(results.Specialization_dprime)

	CIp.Specialization_pred_dprime <- rowQuantiles(results.Specialization_pred_dprime, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, drop = FALSE)	## Na's are excluded first, do not drop singleton dimensions in the result
		rownames(CIp.Specialization_pred_dprime) <- rownames(results.Specialization_pred_dprime)


## Median per row for Specialization dprime and dprime transpose
	med.Specialization_dprime <- t(t(rowMedians(results.Specialization_dprime, na.rm=TRUE)))			## Na's are excluded first
		rownames(med.Specialization_dprime) <- rownames(results.Specialization_dprime)
	med.Specialization_pred_dprime <- t(t(rowMedians(results.Specialization_pred_dprime, na.rm=TRUE)))	## Na's are excluded first
		rownames(med.Specialization_pred_dprime) <- rownames(results.Specialization_pred_dprime)

## Mean per row for Specialization dprime and dprime transpose
	mean.Specialization_dprime <- t(t(rowMeans(results.Specialization_dprime, na.rm=TRUE)))			## Na's are excluded first
		rownames(mean.Specialization_dprime) <- rownames(results.Specialization_dprime)
	mean.Specialization_pred_dprime <- t(t(rowMeans(results.Specialization_pred_dprime, na.rm=TRUE)))	## Na's are excluded first
		rownames(mean.Specialization_pred_dprime) <- rownames(results.Specialization_pred_dprime)

## Contingency Table Results
	if(data_type == "ALL") {
	med.quant_fisher_central_oddsratio <- median(results.quant_fisher_central_oddsratio)
	med.quant_fisher_central_pval <- median(results.quant_fisher_central_pval)
	med.quant_fisher_minlike_oddsratio <- median(results.quant_fisher_minlike_oddsratio)
	med.quant_fisher_minlike_pval <- median(results.quant_fisher_minlike_pval)
	med.quant_blaker_oddsratio <- median(results.quant_blaker_oddsratio)
	med.quant_blaker_pval <- median(results.quant_blaker_pval)

	med.ind_fisher_central_oddsratio <- median(results.ind_fisher_central_oddsratio)
	med.ind_fisher_central_pval <- median(results.ind_fisher_central_pval)
	med.ind_fisher_minlike_oddsratio <- median(results.ind_fisher_minlike_oddsratio)
	med.ind_fisher_minlike_pval <- median(results.ind_fisher_minlike_pval)
	med.ind_blaker_oddsratio <- median(results.ind_blaker_oddsratio)
	med.ind_blaker_pval <- median(results.ind_blaker_pval)

	CIp.quant_fisher_central_oddsratio <- as.numeric(quantile(results.quant_fisher_central_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.quant_fisher_central_pval <- as.numeric(quantile(results.quant_fisher_central_pval, c(alpha / 2, 1 - alpha / 2)))
	CIp.quant_fisher_minlike_oddsratio <- as.numeric(quantile(results.quant_fisher_minlike_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.quant_fisher_minlike_pval <- as.numeric(quantile(results.quant_fisher_minlike_pval, c(alpha / 2, 1 - alpha / 2)))
	CIp.quant_blaker_oddsratio <- as.numeric(quantile(results.quant_blaker_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.quant_blaker_pval <- as.numeric(quantile(results.quant_blaker_pval, c(alpha / 2, 1 - alpha / 2)))

	CIp.ind_fisher_central_oddsratio <- as.numeric(quantile(results.ind_fisher_central_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.ind_fisher_central_pval <- as.numeric(quantile(results.ind_fisher_central_pval, c(alpha / 2, 1 - alpha / 2)))
	CIp.ind_fisher_minlike_oddsratio <- as.numeric(quantile(results.ind_fisher_minlike_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.ind_fisher_minlike_pval <- as.numeric(quantile(results.ind_fisher_minlike_pval, c(alpha / 2, 1 - alpha / 2)))
	CIp.ind_blaker_oddsratio <- as.numeric(quantile(results.ind_blaker_oddsratio, c(alpha / 2, 1 - alpha / 2)))
	CIp.ind_blaker_pval <- as.numeric(quantile(results.ind_blaker_pval, c(alpha / 2, 1 - alpha / 2)))
	} ## END CONTINGENCY TABLE RESULTS
	
################ END OF BOOTSTRAP RESULTS
################
################
################
################
################
################
################
################ ORIGINAL DATASET EXTRA RESULTS
	## Group Spiders Together to Generate Quantitative Data
		for (i in colnames(dataset)){
			sql1 <- fn$identity("
				Select Family_Short,
					sum(Agel) as Agel,
					sum(Cori) as Cori,
					sum(Dict) as Dict,
					sum(Gnap) as Gnap,
					sum(Hahn) as Hahn,
					sum(Liny) as Liny,
					sum(Lyco) as Lyco,
					sum(Pisa) as Pisa,
					sum(Salt) as Salt,
					sum(Ther) as Ther,
					sum(Thom) as Thom,
					sum(Arac) as Arac,
					sum(Cole) as Cole,
					sum(Coll) as Coll,
					sum(Derm) as Derm,
					sum(Dip) as Dip,
					sum(Gryl) as Gryl,
					sum(Hyme) as Hyme,
					sum(Isop) as Isop,
					sum(Lepi) as Lepi,
					sum(Opil) as Opil,
					sum(Pseu) as Pseu,
					count(Family_Short) as Total_Spiders
				From dataset
				Group by Family_Short ")
			dataset_grouped <- sqldf(c(sql1))					## Grouped Dataset
				}

		orig_Quantitative_frame <- dataset_grouped[2:24]			## Quantitative Dataset as Data Frame
		row.names(orig_Quantitative_frame) <- t(dataset_grouped[1])		## Set the Row Names as the Family Names
		orig_Quantitative_matrix <- as.matrix(orig_Quantitative_frame)	## Quantitative Dataset as Matrix


	## Generate Indicator Dataset without X's where RowName = ColName (also removes total number of spiders column)
		orig_Indicator_frame <- orig_Quantitative_frame[1:22]
		orig_Indicator_frame_rownames <- row.names(orig_Indicator_frame)

		## Use SQL to Change All Values > 0 to 1
			for (i in colnames(orig_Indicator_frame)) {
				sql1 <- fn$identity("
					Update orig_Indicator_frame
					Set $i = '$Link.yes'
					Where $i > '$Link.no' ")
				sql2 <- "Select * From main.orig_Indicator_frame"
				orig_Indicator_frame <- sqldf(c(sql1, sql2))		## Indicator Dataset as Data Frame
				}

			## Restore Row Names
			row.names(orig_Indicator_frame) <- orig_Indicator_frame_rownames

			## Create Indicator Dataset as Matrix
			orig_Indicator_matrix <- as.matrix(orig_Indicator_frame)	## Indicator Dataset as Matrix

	## Undirected Indicator Dataset with All Predators (For use with Connectance)
		temp3 <- orig_Indicator_frame

		## Add Predators that Don't currently exist in the Directed Dataset
		if(rownames(temp3) %contains% "Agel" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Agel")
			}

		if(rownames(temp3) %contains% "Cori" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Cori")
			}

		if(rownames(temp3) %contains% "Dict" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Dict")
			}

		if(rownames(temp3) %contains% "Gnap" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Gnap")
			}

		if(rownames(temp3) %contains% "Hahn" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Hahn")
			}

		if(rownames(temp3) %contains% "Liny" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Liny")
			}

		if(rownames(temp3) %contains% "Lyco" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Lyco")
			}

		if(rownames(temp3) %contains% "Pisa" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Pisa")
			}

		if(rownames(temp3) %contains% "Salt" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Salt")
			}

		if(rownames(temp3) %contains% "Ther" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Ther")
			}

		if(rownames(temp3) %contains% "Thom" == FALSE){
			temp3_row_names_temp <- row.names(temp3)
			temp3<- rbind(temp3,vec0)
			row.names(temp3) <- c(temp3_row_names_temp,"Thom")
			}

		## Sort Dataset by Family and also renames rows
		temp3_row_names <- row.names(temp3)			## Generates Final List of Row Names
		temp3["Family_Short"] <- temp3_row_names		## Stores List of Row Names in Family_Short Column

		for (i in colnames(temp3)){
			sql1 <- fn$identity("
				Select * From temp3
				Order by Family_Short asc ")
			temp3_sorted <- sqldf(c(sql1))		## Sorts temp1
			}

		row.names(temp3_sorted) <- t(temp3_sorted["Family_Short"])	## Gives Row Names to Sorted Temp3

		temp4 <- temp3_sorted[1:22]				## Gets rid of Family_Short Column in Temp3


		## Verify that Row Names Match First Predator Column Names
		if(is.logical(rownames(temp4) == colnames(temp4[1:numberpredatoritems])) == FALSE){break}

		## Create Undirected Indicator Matrix (for use with Pimm & Lawton, and Undirected Connectance)
		orig_Indicator_undirected_matrix <- matrix(data = 0, nrow = ncol(temp4), ncol = ncol(temp4))

		## Update Undirected Indicator Matrix so that all links are shown
		## Predators = rows, Prey = columns
		## Predators = k, Prey = l
		for(k in 1:nrow(temp4)){
			for(l in 1:ncol(temp4)){
				if(temp4[k,l] == 1){(orig_Indicator_undirected_matrix[k,l] <- 1) & (orig_Indicator_undirected_matrix[l,k] <- 1)}
			}
		}

		## Test to verify if symmetric or not.  If not symmetric, then process stops because of error.
			if(isSymmetric(orig_Indicator_undirected_matrix) == FALSE){break}

		## Add Back Row Names
		rownames(orig_Indicator_undirected_matrix) <- colnames(temp4)
		colnames(orig_Indicator_undirected_matrix) <- colnames(temp4)


	## Undirected Connectance
		orig_L <- sum(colSums(orig_Indicator_undirected_matrix))
		orig_S <- numberpredatoritems + numberpreyitems
		orig_Connectance_wcannibalism <- orig_L / orig_S^2
		orig_Connectance_undirected <- orig_L / ( orig_S * (orig_S - 1) )

	## IGP Undirected Connectance Predators Only
		orig_IGP_Conn_Pred_undirected <- orig_Indicator_undirected_matrix[1:numberpredatoritems, 1:numberpredatoritems]
		orig_L_IGP_pred <- sum(colSums(orig_IGP_Conn_Pred_undirected))
		orig_S_IGP_pred <- numberpredatoritems
		orig_Connectance_undirected_pred <- orig_L_IGP_pred / ( orig_S_IGP_pred * (orig_S_IGP_pred - 1) )


	## Contingency Table
	if(data_type == "ALL"){
		## Splitting Datasets
			for(i in colnames(dataset_grouped)){
				sql1 <- fn$identity("
					Select Family_Short,
						sum(Agel) as Agel,
						sum(Cori) as Cori,
						sum(Dict) as Dict,
						sum(Gnap) as Gnap,
						sum(Hahn) as Hahn,
						sum(Liny) as Liny,
						sum(Lyco) as Lyco,
						sum(Pisa) as Pisa,
						sum(Salt) as Salt,
						sum(Ther) as Ther,
						sum(Thom) as Thom
					From dataset_grouped
					Where Family_Short in ('Cori', 'Gnap', 'Lyco', 'Pisa', 'Salt', 'Thom')
					Group by Family_Short ")
				sql2 <- fn$identity("
					Select Family_Short,
						sum(Arac) as Arac,
						sum(Cole) as Cole,
						sum(Coll) as Coll,
						sum(Derm) as Derm,
						sum(Dip) as Dip,
						sum(Gryl) as Gryl,
						sum(Hyme) as Hyme,
						sum(Isop) as Isop,
						sum(Lepi) as Lepi,
						sum(Opil) as Opil,
						sum(Pseu) as Pseu
					From dataset_grouped
					Where Family_Short in ('Cori', 'Gnap', 'Lyco', 'Pisa', 'Salt', 'Thom')
					Group by Family_Short ")
				sql3 <- fn$identity("
					Select Family_Short,
						sum(Agel) as Agel,
						sum(Cori) as Cori,
						sum(Dict) as Dict,
						sum(Gnap) as Gnap,
						sum(Hahn) as Hahn,
						sum(Liny) as Liny,
						sum(Lyco) as Lyco,
						sum(Pisa) as Pisa,
						sum(Salt) as Salt,
						sum(Ther) as Ther,
						sum(Thom) as Thom
					From dataset_grouped
					Where Family_Short in ('Agel', 'Dict', 'Hahn', 'Liny', 'Ther')
					Group by Family_Short ")
				sql4 <- fn$identity("
					Select Family_Short,
						sum(Arac) as Arac,
						sum(Cole) as Cole,
						sum(Coll) as Coll,
						sum(Derm) as Derm,
						sum(Dip) as Dip,
						sum(Gryl) as Gryl,
						sum(Hyme) as Hyme,
						sum(Isop) as Isop,
						sum(Lepi) as Lepi,
						sum(Opil) as Opil,
						sum(Pseu) as Pseu
					From dataset_grouped
					Where Family_Short in ('Agel', 'Dict', 'Hahn', 'Liny', 'Ther')
					Group by Family_Short ")
				orig_dataset_curs_spiders <- sqldf(c(sql1))
				orig_dataset_curs_nonspiders <- sqldf(c(sql2))
				orig_dataset_web_spiders <- sqldf(c(sql3))
				orig_dataset_web_nonspiders <- sqldf(c(sql4))
				}


				orig_dataset_curs_spiders2 <- orig_dataset_curs_spiders[,2:(numberpredatoritems+1)]
				orig_dataset_curs_nonspiders2 <- orig_dataset_curs_nonspiders[,2:(numberpreyitems+1)]
				orig_dataset_web_spiders2 <- orig_dataset_web_spiders[,2:(numberpredatoritems+1)]
				orig_dataset_web_nonspiders2 <- orig_dataset_web_nonspiders[,2:(numberpreyitems+1)]

				row.names(orig_dataset_curs_spiders2) <- t(orig_dataset_curs_spiders["Family_Short"])
				row.names(orig_dataset_curs_nonspiders2) <- t(orig_dataset_curs_nonspiders["Family_Short"])
				row.names(orig_dataset_web_spiders2) <- t(orig_dataset_web_spiders["Family_Short"])
				row.names(orig_dataset_web_nonspiders2) <- t(orig_dataset_web_nonspiders["Family_Short"])


		## Changing Split Datasets to Indicators
			orig_indicator_curs_spiders <- orig_dataset_curs_spiders2
			orig_indicator_curs_nonspiders <- orig_dataset_curs_nonspiders2
			orig_indicator_web_spiders <- orig_dataset_web_spiders2
			orig_indicator_web_nonspiders <- orig_dataset_web_nonspiders2

			orig_indicator_curs_spiders_rownames <- row.names(orig_indicator_curs_spiders)
			orig_indicator_curs_nonspiders_rownames <- row.names(orig_indicator_curs_nonspiders)
			orig_indicator_web_spiders_rownames <- row.names(orig_indicator_web_spiders)
			orig_indicator_web_nonspiders_rownames <- row.names(orig_indicator_web_nonspiders)


			for (i in colnames(orig_indicator_curs_spiders)) {
				sql1 <- fn$identity("
					Update orig_indicator_curs_spiders
					Set $i = '$Link.yes'
					Where $i > '$Link.no' ")
				sql2 <- "Select * From orig_indicator_curs_spiders"
				orig_indicator_curs_spiders <- sqldf(c(sql1,sql2))
				}

			for (i in colnames(orig_indicator_curs_nonspiders)) {
				sql1 <- fn$identity("
					Update orig_indicator_curs_nonspiders
					Set $i = '$Link.yes'
					Where $i > '$Link.no' ")
				sql2 <- "Select * From orig_indicator_curs_nonspiders"
				orig_indicator_curs_nonspiders <- sqldf(c(sql1,sql2))
				}

			for (i in colnames(orig_indicator_web_spiders)) {
				sql1 <- fn$identity("
					Update orig_indicator_web_spiders
					Set $i = '$Link.yes'
					Where $i > '$Link.no' ")
				sql2 <- "Select * From orig_indicator_web_spiders"
				orig_indicator_web_spiders <- sqldf(c(sql1,sql2))
				}

			for (i in colnames(orig_indicator_web_nonspiders)) {
				sql1 <- fn$identity("
					Update orig_indicator_web_nonspiders
					Set $i = '$Link.yes'
					Where $i > '$Link.no' ")
				sql2 <- "Select * From orig_indicator_web_nonspiders"
				orig_indicator_web_nonspiders <- sqldf(c(sql1,sql2))
				}

			row.names(orig_indicator_curs_spiders) <- orig_indicator_curs_spiders_rownames
			row.names(orig_indicator_curs_nonspiders) <- orig_indicator_curs_nonspiders_rownames
			row.names(orig_indicator_web_spiders) <- orig_indicator_web_spiders_rownames
			row.names(orig_indicator_web_nonspiders) <- orig_indicator_web_nonspiders_rownames


		## Count of Spiders Quant
			orig_quant_count_curs_spiders <- sum(orig_dataset_curs_spiders2)
			orig_quant_count_curs_nonspiders <- sum(orig_dataset_curs_nonspiders2)
			orig_quant_count_web_spiders <- sum(orig_dataset_web_spiders2)
			orig_quant_count_web_nonspiders <- sum(orig_dataset_web_nonspiders2)


		## Count of Spiders Indicator
			orig_ind_count_curs_spiders <- sum(orig_indicator_curs_spiders)
			orig_ind_count_curs_nonspiders <- sum(orig_indicator_curs_nonspiders)
			orig_ind_count_web_spiders <- sum(orig_indicator_web_spiders)
			orig_ind_count_web_nonspiders <- sum(orig_indicator_web_nonspiders)

		## Contingency Table Quant
			orig_quant_cont_table <- matrix(c(
					orig_quant_count_curs_spiders, orig_quant_count_curs_nonspiders,
					orig_quant_count_web_spiders, orig_quant_count_web_nonspiders),
				nrow = 2, ncol = 2,
				byrow = TRUE,
				dimnames = list(c("Cursorial","Web Building"),
							c("Spiders", "Non-Spiders")))

			orig_quant_fisher_central <- fisher.exact(orig_quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "central")
			orig_quant_fisher_minlike <- fisher.exact(orig_quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "minlike")
			orig_quant_blaker <- blaker.exact(orig_quant_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05)


		## Contingency Table Indicator
			orig_ind_cont_table <- matrix(c(
					orig_ind_count_curs_spiders, orig_ind_count_curs_nonspiders,
					orig_ind_count_web_spiders, orig_ind_count_web_nonspiders),
				nrow = 2, ncol = 2,
				byrow = TRUE,
				dimnames = list(c("Cursorial", "Web Building"),
							c("Spiders", "Non-Spiders")))

			orig_ind_fisher_central <- fisher.exact(orig_ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "central")
			orig_ind_fisher_minlike <- fisher.exact(orig_ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05, tsmethod = "minlike")
			orig_ind_blaker <- blaker.exact(orig_ind_cont_table, alternative = "two.sided", conf.int = TRUE, conf.level = 1-alpha, tol = 1e-05)
		} ## END OF CONTINGENCY TABLE

	## Cleanup
		rm(sql1)
		rm(sql2)
		rm(sql3)
		rm(sql4)
		rm(sql5)
		rm(sql6)
		rm(sql7)
		rm(sql8)
		rm(sql9)
		rm(sql10)
		rm(sql11)
		rm(temp3)
		rm(temp3_row_names)
		rm(temp3_sorted)

################ END OF ORIGINAL DATASET CALCULATIONS
################
################
################
################
################
################ Print Results
	print("Percentile Bootstrap Confidence Intervals")
	print(CIp.Connectance_wcannibalism)		## Undirected Connectance with cannibalism (denominator S^2)
	print(CIp.Connectance_undirected)		## Undirected Connectance
	print(CIp.Connectance_directed)		## Directed Connectance
	print(CIp.Interaction_Strength)		## Interaction Strength
	print(CIp.IGP_Conn_Pred_undirected)		## IGP Undirected Connectance Predators
	print(CIp.IGP_Conn_Pred_directed)		## IGP Directed Connectance Predators
	print(CIp.IGP_Conn_Prey_directed)		## IGP Directed Connectance Prey
	print(CIp.IGP_IS_Pred)				## IGP Interaction Strength Predators
	print(CIp.IGP_IS_Prey)				## IGP Interaction Strength Prey
	print(CIp.Interaction_Evenness_e_sum)	## Interaction Evenness; logbase = e; intereven = sum
	print(CIp.Interaction_Evenness_2_sum)	## Interaction Evenness; logbase = 2; intereven = sum
	print(CIp.Interaction_Evenness_e_prod)	## Interaction Evenness; logbase = e; intereven = prod
	print(CIp.Interaction_Evenness_2_prod)	## Interaction Evenness; logbase = 2; intereven = prod
	print(CIp.Specialization_dprime)		## Specialization dprime; Pred = col, Prey = row
	print(CIp.Specialization_pred_dprime)	## Specialization dprime transpose; Prey = col, Pred = row
	print(CIp.Specialization_h2)			## Specialization h2; Pred = col, Prey = row
	print(CIp.Specialization_h2_uncorr)		## Specialization h2 uncorrelated; Pred = col, Prey = row
	print(CIp.Specialization_pred_h2)		## Specialization h2; Prey = col, Pred = row
	print(CIp.Specialization_pred_h2_uncorr)	## Specialization h2 uncorrelated; Prey = col, Pred = row
	print(CIp.Compart_Bipartite)			## Compartmentalization using Bipartite
	print(CIp.Compartmentalization_pimm_all)	## Pimm & Lawton Compartmentalization; All Spiders
	print(CIp.Compartmentalization_pimm_pred)	## Pimm & Lawton Compartmentalization; Predators Only
	


	print("Median Values")
	print(med.Connectance_wcannibalism)		## Undirected Connectance with cannibalism (denominator S^2)
	print(med.Connectance_undirected)		## Undirected Connectance
	print(med.Connectance_directed)		## Directed Connectance
	print(med.Interaction_Strength)		## Interaction Strength
	print(med.IGP_Conn_Pred_undirected)		## IGP Undirected Connectance Predators
	print(med.IGP_Conn_Pred_directed)		## IGP Directed Connectance Predators
	print(med.IGP_Conn_Prey_directed)		## IGP Directed Connectance Prey
	print(med.IGP_IS_Pred)				## IGP Interaction Strength Predators
	print(med.IGP_IS_Prey)				## IGP Interaction Strength Prey
	print(med.Interaction_Evenness_e_sum)	## Interaction Evenness; logbase = e; intereven = sum
	print(med.Interaction_Evenness_2_sum)	## Interaction Evenness; logbase = 2; intereven = sum
	print(med.Interaction_Evenness_e_prod)	## Interaction Evenness; logbase = e; intereven = prod
	print(med.Interaction_Evenness_2_prod)	## Interaction Evenness; logbase = 2; intereven = prod
	print(med.Specialization_dprime)		## Specialization dprime; Pred = col, Prey = row
	print(med.Specialization_pred_dprime)	## Specialization dprime transpose; Prey = col, Pred = row
	print(med.Specialization_h2)			## Specialization h2; Pred = col, Prey = row
	print(med.Specialization_h2_uncorr)		## Specialization h2 uncorrelated; Pred = col, Prey = row
	print(med.Specialization_pred_h2)		## Specialization h2; Prey = col, Pred = row
	print(med.Specialization_pred_h2_uncorr)	## Specialization h2 uncorrelated; Prey = col, Pred = row
	print(med.Compart_Bipartite)			## Compartmentalization using Bipartite
	print(med.Compartmentalization_pimm_all)	## Pimm & Lawton Compartmentalization; All Spiders
	print(med.Compartmentalization_pimm_pred)	## Pimm & Lawton Compartmentalization; Predators Only



	if(data_type == "ALL"){
		print("Contingency Table Results Quant Data - Median Odds Ratios")
		print(med.quant_fisher_central_oddsratio)
		print(med.quant_fisher_minlike_oddsratio)
		print(med.quant_blaker_oddsratio)


		print("Contingency Table Results Quant Data - Median P Values")
		print(med.quant_fisher_central_pval)
		print(med.quant_fisher_minlike_pval)
		print(med.quant_blaker_pval)


		print("Contingency Table Results Indicator Data - Median Odds Ratios")
		print(med.ind_fisher_central_oddsratio)
		print(med.ind_fisher_minlike_oddsratio)
		print(med.ind_blaker_oddsratio)


		print("Contingency Table Results Indicator Data - Median P Values")
		print(med.ind_fisher_central_pval)
		print(med.ind_fisher_minlike_pval)
		print(med.ind_blaker_pval)



		print("Contingency Table Results Quant Data - Confidence Interval on the Odds Ratios")
		print(CIp.quant_fisher_central_oddsratio)
		print(CIp.quant_fisher_minlike_oddsratio)
		print(CIp.quant_blaker_oddsratio)



		print("Contingency Table Results Quant Data - Confidence Interval on the P Values")
		print(CIp.quant_fisher_central_pval)
		print(CIp.quant_fisher_minlike_pval)
		print(CIp.quant_blaker_pval)



		print("Contingency Table Results Indicator Data - Confidence Interval on the Odds Ratios")
		print(CIp.ind_fisher_central_oddsratio)
		print(CIp.ind_fisher_minlike_oddsratio)
		print(CIp.ind_blaker_oddsratio)



		print("Contingency Table Results Indicator Data - Confidence Interval on the P Values")
		print(CIp.ind_fisher_central_pval)
		print(CIp.ind_fisher_minlike_pval)
		print(CIp.ind_blaker_pval)
	} ## END PRINT WHERE DATA_TYPE = ALL


	print("Original Dataset Extra Results")
	print(c(orig_Connectance_wcannibalism, orig_L, orig_S^2))							## Undirected Connectance with cannibalism, Numerator, Denominator
	print(c(orig_Connectance_undirected, orig_L, orig_S * (orig_S - 1)))						## Undirected Connectance, Numerator, Denominator
	print(c(orig_Connectance_undirected_pred, orig_L_IGP_pred, orig_S_IGP_pred * (orig_S_IGP_pred - 1)))	## IGP Undirected Connectance Predators, Numerator, Denominator


	if(data_type == "ALL") {
		print("Original Dataset Contingency Table Indicator Results")
		print(orig_ind_fisher_central)		## Indicator Data; Central Fisher Exact Test
		print(orig_ind_fisher_minlike)		## Indicator Data; Fisher Exact Test with Minimum Likelihood
		print(orig_ind_blaker)				## Indicator Data; Blaker Exact Test


		print("Original Dataset Contingency Table Quant Results")
		print(orig_quant_fisher_central)		## Quant Data; Central Fisher Exact Test
		print(orig_quant_fisher_minlike)		## Quant Data; Fisher Exact Test with Minimum Likelihood
		print(orig_quant_blaker)			## Quant Data; Blaker Exact Test
	}


## Sample Means
	mean.Connectance_wcannibalism <- mean(results.Connectance_wcannibalism)
	mean.Connectance_undirected <- mean(results.Connectance_undirected)
	mean.Connectance_directed <- mean(results.Connectance_directed)
	mean.Interaction_Strength <- mean(results.Interaction_Strength)
	mean.IGP_Conn_Pred_undirected <- mean(results.IGP_Conn_Pred_undirected)
	mean.IGP_Conn_Pred_directed <- mean(results.IGP_Conn_Pred_directed)
	mean.IGP_Conn_Prey_directed <- mean(results.IGP_Conn_Prey_directed)
	mean.IGP_IS_Pred <- mean(results.IGP_IS_Pred)
	mean.IGP_IS_Prey <- mean(results.IGP_IS_Prey)
	mean.Interaction_Evenness_e_sum <- mean(results.Interaction_Evenness_e_sum)
	mean.Interaction_Evenness_2_sum <- mean(results.Interaction_Evenness_2_sum)
	mean.Interaction_Evenness_e_prod <- mean(results.Interaction_Evenness_e_prod)
	mean.Interaction_Evenness_2_prod <- mean(results.Interaction_Evenness_2_prod)
	mean.Specialization_h2 <- mean(results.Specialization_h2)
	mean.Specialization_h2_uncorr <- mean(results.Specialization_h2_uncorr)
	mean.Specialization_pred_h2 <- mean(results.Specialization_pred_h2)
	mean.Specialization_pred_h2_uncorr <- mean(results.Specialization_pred_h2_uncorr)
	mean.Compart_Bipartite <- mean(results.Compart_Bipartite)
	mean.Compartmentalization_pimm_all <- mean(results.Compartmentalization_pimm_all)
	mean.Compartmentalization_pimm_pred <- mean(results.Compartmentalization_pimm_pred)

	print("Mean Values")
	print(mean.Connectance_wcannibalism)		## Undirected Connectance with cannibalism (denominator S^2)
	print(mean.Connectance_undirected)		## Undirected Connectance
	print(mean.Connectance_directed)		## Directed Connectance
	print(mean.Interaction_Strength)		## Interaction Strength
	print(mean.IGP_Conn_Pred_undirected)		## IGP Undirected Connectance Predators
	print(mean.IGP_Conn_Pred_directed)		## IGP Directed Connectance Predators
	print(mean.IGP_Conn_Prey_directed)		## IGP Directed Connectance Prey
	print(mean.IGP_IS_Pred)				## IGP Interaction Strength Predators
	print(mean.IGP_IS_Prey)				## IGP Interaction Strength Prey
	print(mean.Interaction_Evenness_e_sum)	## Interaction Evenness; logbase = e; intereven = sum
	print(mean.Interaction_Evenness_2_sum)	## Interaction Evenness; logbase = 2; intereven = sum
	print(mean.Interaction_Evenness_e_prod)	## Interaction Evenness; logbase = e; intereven = prod
	print(mean.Interaction_Evenness_2_prod)	## Interaction Evenness; logbase = 2; intereven = prod
	print(mean.Specialization_dprime)		## Specialization dprime; Pred = col, Prey = row
	print(mean.Specialization_pred_dprime)	## Specialization dprime transpose; Prey = col, Pred = row
	print(mean.Specialization_h2)			## Specialization h2; Pred = col, Prey = row
	print(mean.Specialization_h2_uncorr)		## Specialization h2 uncorrelated; Pred = col, Prey = row
	print(mean.Specialization_pred_h2)		## Specialization h2; Prey = col, Pred = row
	print(mean.Specialization_pred_h2_uncorr)	## Specialization h2 uncorrelated; Prey = col, Pred = row
	print(mean.Compart_Bipartite)			## Compartmentalization using Bipartite
	print(mean.Compartmentalization_pimm_all)	## Pimm & Lawton Compartmentalization; All Spiders
	print(mean.Compartmentalization_pimm_pred)	## Pimm & Lawton Compartmentalization; Predators Only





