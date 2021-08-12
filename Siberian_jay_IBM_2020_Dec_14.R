####################################################################################################################################################################
#################################################################### IBM Siberian jay ##############################################################################
########################################################## Written by Ute Bradter, ute.bradter@nina.no  ############################################################
############################################################ Version: FunIBM6, 14 Dec 2020 #########################################################################				
####################################################################################################################################################################

library(sp)
library(rgdal)
library(raster)
library(rgeos)

################################################################################### IBM function ###################################################################
####################################################################################################################################################################

FunSJIBM <- function(YearMin, YearMax, YearsScenarioPeriod, StandsShp = NULL, CovarDir, StartNoRep, NoRep, SJLogInit, TargetAreaName, NamingComment, TempFolder, OutFolder){

	for (Rep in (StartNoRep : NoRep)){

# Simulation years
		Years <- seq(YearMin, YearMax, 1)

# Scenario years (a new forest scenario created at YearsScenarioPeriod intervals
		ScenarioYears <- rep(c((YearMin - 1) : round(length(Years)/YearsScenarioPeriod-1)), each = YearsScenarioPeriod)
		
# create one file of survival and transition betas per stage for predicting via model matrix
		MarkBetas <- BestNbh[[1]]$results$beta

		BetasE <- MarkBetas[c(1, 3, 7, 15), 1]
		names(BetasE) <- rownames(MarkBetas[c(1, 3, 7, 15),])
		names(BetasE) <- gsub("S:|:B3", "", names(BetasE))
	
		BetasF <- MarkBetas[c(1, 4, 7, 10, 14, 18), 1]
		names(BetasF) <- rownames(MarkBetas[c(1, 4, 7, 10, 14, 18),])
		names(BetasF) <- gsub("S:|S:B9:|:B9", "", names(BetasF))

		BetasC <- MarkBetas[c(1, 7), 1]
		names(BetasC) <- rownames(MarkBetas[c(1, 7),])
		names(BetasC) <- gsub("S:", "", names(BetasC))

		BetasD <- MarkBetas[c(1, 2, 7, 9, 13, 17), 1]
		names(BetasD) <- rownames(MarkBetas[c(1, 2, 7, 9, 13, 17),])
		names(BetasD) <- gsub("S:|S:NB9:|:NB9", "", names(BetasD))

		BetasX <- MarkBetas[c(1, 5, 7), 1]
		names(BetasX) <- rownames(MarkBetas[c(1, 5, 7),])
		names(BetasX) <- gsub("S:", "", names(BetasX))

		BetasZ <- MarkBetas[c(1, 6, 7, 8, 11, 12, 16, 19), 1]
		names(BetasZ) <- rownames(MarkBetas[c(1, 6, 7, 8, 11, 12, 16, 19),])
		names(BetasZ) <- gsub("S:|S:DisJ:|:DisJ", "", names(BetasZ))

		BetasNB3ToB9 <- MarkBetas[c(25, 29, 30, 39), 1]
		names(BetasNB3ToB9) <- rownames(MarkBetas[c(25, 29, 30, 39),])
		names(BetasNB3ToB9) <- gsub("Psi:|:NB3ToB9", "", names(BetasNB3ToB9))

		BetasNB9ToB3 <- MarkBetas[c(25, 28, 29, 30, 33, 38, 40), 1]
		names(BetasNB9ToB3) <- rownames(MarkBetas[c(25, 28, 29, 30, 33, 38, 40),])
		names(BetasNB9ToB3) <- gsub("Psi:", "", names(BetasNB9ToB3))
		names(BetasNB9ToB3) <- gsub("NB9ToB3:", "", names(BetasNB9ToB3))

		BetasRJToB3 <- MarkBetas[c(25, 26, 29, 30, 31, 34, 36), 1]
		names(BetasRJToB3) <- rownames(MarkBetas[c(25, 26, 29, 30, 31, 34, 36),])
		names(BetasRJToB3) <- gsub("Psi:", "", names(BetasRJToB3))
		names(BetasRJToB3) <- gsub("RJToB3:", "", names(BetasRJToB3))

		BetasDJToB3 <- MarkBetas[c(25, 27, 29, 30, 32, 35, 37), 1]
		names(BetasDJToB3) <- rownames(MarkBetas[c(25, 27, 29, 30, 32, 35, 37),])
		names(BetasDJToB3) <- gsub("Psi:", "", names(BetasDJToB3))
		names(BetasDJToB3) <- gsub("DJToB3:", "", names(BetasDJToB3))

# a file to log Mark mean estimates
		MarkMeanEstimate <- as.data.frame(matrix(nrow = length(Years), ncol = 11))
		colnames(MarkMeanEstimate) <- c("year", "S_E", "S_C", "S_F", "S_D", "S_X", "S_Z", "Psi_CF", "Psi_DE", "Psi_XE", "Psi_ZE")	
		MarkMeanEstimate$year <- Years
# a file to keep a tally
		SJTally <- as.data.frame(matrix(ncol = 16, nrow = (length(Years) * 4 - 2)))
		colnames(SJTally) <- c("Year", "Season", "Process", "RetJuv", "DispJuv", "B", "NB", "TransNB", "TransRJ", "TransDJ", "EmiB", "EmiNB", "MeanGSTerr", "MaxGSTerr", "NoTerrs", "Only1B")
		SJTally$Year <- c(rep(Years[1], each = 2), rep(Years[2 : length(Years)], each = 4))
		SJTally$Season <- c(rep(c("Sept", "Sept", "March", "March"), (length(Years)-1)), "Sept", "Sept")
		SJTally$Process <- rep(c("AfterSurv", "AfterAll"), nrow(SJTally) / 2)
# start of the log file
		SJLog <- SJLogInit

		SJTallyInit <- as.data.frame(matrix(ncol = 16, nrow = 1))
		colnames(SJTallyInit) <- c("Year", "Season", "Process", "RetJuv", "DispJuv", "B", "NB", "TransNB", "TransRJ", "TransDJ", "EmiB", "EmiNB", "MeanGSTerr", "MaxGSTerr", "NoTerrs", "Only1B")
		SJTallyInit$Year <- YearMin
		SJTallyInit$Season <- "March"
		SJTallyInit$Process <- "AfterAll"
		SJTallyInit$B <- sum(SJLogInit$E)
		SJTallyInit$NB <- sum(SJLogInit$C)
		SJTallyInit$MeanGSTerr <- mean(apply(SJLogInit[, c("E", "F", "C", "D", "X", "Z")], 1, sum))
		SJTallyInit$MaxGSTerr <- max(apply(SJLogInit[, c("E", "F", "C", "D", "X", "Z")], 1, sum))
		SJTallyInit$NoTerrs <- nrow(SJLogInit)
		SJTallyInit$Only1B <- nrow(SJLogInit[SJLogInit$E == 1, ])

		for (y in (1 : length(Years))){
			
# covars at 30 ha scale on a 250 m grid
			Preds30ha <- read.csv(paste(CovarDir, "\\PercSpruce30ha", NamingComment, ScenarioYears[y], ".csv", sep = ""))
			colnames(Preds30ha)[which(colnames(Preds30ha) == "ID")] <- "GridID"
# covars at 1 km scale on a 250 m grid
			Preds1km <- read.csv(paste(CovarDir, "\\Preds1km", NamingComment, ScenarioYears[y], ".csv", sep = ""))
			colnames(Preds1km)[which(colnames(Preds1km) == "ID")] <- "GridID"
# covars at 2 km scale on a 250 m grid
			Preds2km <- read.csv(paste(CovarDir, "\\Preds2kmIBM", NamingComment, ScenarioYears[y], ".csv", sep = ""))
			colnames(Preds2km)[which(colnames(Preds2km) == "ID")] <- "GridID"
# covars at 10 km scale on a 250 m grid
			Preds10km <- read.csv(paste(CovarDir, "\\Preds10km", NamingComment, ScenarioYears[y], ".csv", sep = ""))
			colnames(Preds10km)[which(colnames(Preds10km) == "ID")] <- "GridID"
			colnames(Preds10km)[which(colnames(Preds10km) == "Nbh")] <- "PatchVAMeanAge10km"

# only points which are in the stands shapefile
			if (length(StandsShp) > 0){	
				MyCoords <- cbind(Preds30ha$x, Preds30ha$y)
				PredsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(Preds30ha), proj4string = CRS("+init=epsg:3021"))
				PredsSPDF <- PredsSPDF[StandsShp, ]
				Preds30ha <- data.frame(PredsSPDF)

				MyCoords <- cbind(Preds1km$x, Preds1km$y)
				PredsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(Preds1km), proj4string = CRS("+init=epsg:3021"))
				PredsSPDF <- PredsSPDF[StandsShp, ]
				Preds1km <- data.frame(PredsSPDF)

				MyCoords <- cbind(Preds2km$x, Preds2km$y)
				PredsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(Preds2km), proj4string = CRS("+init=epsg:3021"))
				PredsSPDF <- PredsSPDF[StandsShp, ]
				Preds2km <- data.frame(PredsSPDF)

				MyCoords <- cbind(Preds10km$x, Preds10km$y)
				PredsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(Preds10km), proj4string = CRS("+init=epsg:3021"))
				PredsSPDF <- PredsSPDF[StandsShp, ]
				Preds10km <- data.frame(PredsSPDF)

				MyCoords <- cbind(PredsStatic$x, PredsStatic$y)
				PredsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(PredsStatic), proj4string = CRS("+init=epsg:3021"))
				PredsSPDF <- PredsSPDF[StandsShp, ]
				PredsStatic <- data.frame(PredsSPDF)
			}

# sort by GridID
			Preds30ha <- Preds30ha[order(Preds30ha$GridID), ]
			Preds1km <- Preds1km[order(Preds1km$GridID), ]
			Preds2km <- Preds2km[order(Preds2km$GridID), ]
			Preds10km <- Preds10km[order(Preds10km$GridID), ]
			PredsStatic <- PredsStatic[order(PredsStatic$GridID), ]

# cap age 
# for demographic model components
			if(max(Preds1km$MAge) > 131){Preds1km[Preds1km$MAge > 131, ]$MAge <- 131}									
			if(max(Preds10km$PatchVAMeanAge10km) > 113){Preds10km[Preds10km$PatchVAMeanAge10km > 113, ]$PatchVAMeanAge10km <- 113}	

# for SDM model component
			if(max(Preds2km$MeanAge) > 169){Preds2km[Preds2km$MeanAge > 169, ]$MeanAge <- 169}	
			
# add covars per model and year
			CurCovarsE <- data.frame(stratumE = 1, NbhAge = Preds10km$PatchVAMeanAge10km, GridID = Preds1km$GridID)
			CurCovarsF <- data.frame(stratumF = 1, NbhAge = Preds10km$PatchVAMeanAge10km, PD = Weather$MeanDaysAboveZero[y], GridID = Preds1km$GridID)
			CurCovarsC <- data.frame(NbhAge = Preds10km$PatchVAMeanAge10km, GridID = Preds1km$GridID)
			CurCovarsD <- data.frame(stratumD = 1, NbhAge = Preds10km$PatchVAMeanAge10km, PD = Weather$MeanDaysAboveZero[y], GridID = Preds1km$GridID)
			CurCovarsX <- data.frame(stratumX = 1, NbhAge = Preds10km$PatchVAMeanAge10km, GridID = Preds1km$GridID)
			CurCovarsZ <- data.frame(stratumZ = 1, NbhAge = Preds10km$PatchVAMeanAge10km, PD = Weather$MeanDaysAboveZero[y], WP = Weather$WinterPrec[y], GridID = Preds1km$GridID)	

			CurCovarsNB3ToB9 <- data.frame(MVol = Preds1km$MVol, GridID = Preds1km$GridID)
			CurCovarsNB9ToB3 <- data.frame(NB9ToB3 = 1, MVol = Preds1km$MVol, PD = Weather$MeanDaysAboveZero[y], GridID = Preds1km$GridID)
			CurCovarsRJToB3 <- data.frame(RJToB3 = 1, MVol = Preds1km$MVol, PD = Weather$MeanDaysAboveZero[y], WP = Weather$WinterPrec[y], GridID = Preds1km$GridID)
			CurCovarsDJToB3 <- data.frame(DJToB3 = 1, MVol = Preds1km$MVol, PD = Weather$MeanDaysAboveZero[y], WP = Weather$WinterPrec[y], GridID = Preds1km$GridID)

			CurCovarsRepro <- data.frame(PatchVAMeanAge10km = Preds10km$PatchVAMeanAge10km, PercSpruce = Preds30ha$PercSpruce, MeanDaysAboveZero = Weather$MeanDaysAboveZero[y], 
				GridID = Preds30ha$GridID)
			CurCovarsEmiBToBMarch <- data.frame(stg3F = 0, PatchVAMeanAge10km = Preds10km$PatchVAMeanAge10km, SpringPrec = Weather$SpringPrec[y], WinterPrec = 0, WinterTemp = 0, 
				GridID = Preds10km$GridID)
			CurCovarsEmiBToBSept <- data.frame(stg3F = 1, PatchVAMeanAge10km = Preds10km$PatchVAMeanAge10km, SpringPrec = 0, WinterPrec = Weather$WinterPrec[y], 
				WinterTemp = Weather$WinterTemp[y], GridID = Preds10km$GridID)
			CurCovarsHS <- data.frame(x = Preds2km$x, y = Preds2km$y, PercVAOF = Preds2km$PercMature, MeanAge = Preds2km$MeanAge, PercVASYF = Preds2km$PercOther, TempJF = PredsStatic$TempJF, 
				PrecAM = PredsStatic$PrecAM, DTM = PredsStatic$DTM, GridID = Preds2km$GridID)

# standardize
			CurCovarsRepro$PatchVAMeanAge10km <- (CurCovarsRepro$PatchVAMeanAge10km - MRepro_PatchVAMeanAge10km) / SRepro_PatchVAMeanAge10km 
			CurCovarsRepro$PercSpruce <- (CurCovarsRepro$PercSpruce - MRepro_PercSpruce) / SRepro_PercSpruce
			CurCovarsRepro$MeanDaysAboveZero <- (CurCovarsRepro$MeanDaysAboveZero - MRepro_PlusDays) / SRepro_PlusDays

			CurCovarsEmiBToBMarch$SpringPrec <- (CurCovarsEmiBToBMarch$SpringPrec - MEmiBToB_Prec) / SEmiBToB_Prec
			CurCovarsEmiBToBSept$WinterPrec <- (CurCovarsEmiBToBSept$WinterPrec - MEmiBToB_Prec) / SEmiBToB_Prec
			CurCovarsEmiBToBSept$WinterTemp <- (CurCovarsEmiBToBSept$WinterTemp - MEmiBToB_Temp) / SEmiBToB_Temp

			CurCovarsHS$PercVAOF <- (CurCovarsHS$PercVAOF - MPercVAOF) /SPercVAOF
			CurCovarsHS$MeanAge <- (CurCovarsHS$MeanAge - MMeanAge) / SMeanAge
			CurCovarsHS$PercVASYF <- (CurCovarsHS$PercVASYF - MPercVASYF) / SPercVASYF
			CurCovarsHS$TempJF <- (CurCovarsHS$TempJF - MTempJF) / STempJF
			CurCovarsHS$PrecAM <- (CurCovarsHS$PrecAM - MPrecAM) / SPrecAM
			CurCovarsHS$DTM <- (CurCovarsHS$DTM - MDTM) / SDTM

##############################################################################################################################################################################
############################################################################ March to September ##############################################################################
##############################################################################################################################################################################

########## the data in March

			if(y == 1){CurSJ <- SJLog} else {CurSJ <- SJLogSpr}
			if(nrow(CurSJ) == 0)break
						
# calculate GSTerr and GSNbh
	
			GSTerrAdd <- FunCalcGSTerr(CurSJ)
			GSNbhAdd <- FunCalcGSNbh(CurSJ)

########## Set up empty log

			CurSJ <- CurSJ[order(CurSJ$GridID), ]
			SJLogAut <- CurSJ
			SJLogAut$Season <- 9
			SJLogAut$Session <- SJLogAut$Session + 1
			SJLogAut[, c("E", "F", "C", "D", "X", "Z")] <- 0

########## 1) Repro model for retained juveniles

# if there are 2 breeders
			CurB <- CurSJ[CurSJ$E > 1, ]
# the covars
			CurCovarsRepro <- CurCovarsRepro[CurCovarsRepro$GridID %in% CurB$GridID, ]
# add group size and standardize
			CurCovarsRepro <- merge(CurCovarsRepro, GSTerrAdd, by = "GridID")
			CurCovarsRepro$GroupSize <- (CurCovarsRepro$GroupSize - MRepro_GroupSize) / SRepro_GroupSize
# order by GridID
			CurCovarsRepro <- CurCovarsRepro[order(CurCovarsRepro$GridID), ]
# simulate from model
			ReproData <- model.matrix(~ PatchVAMeanAge10km + PercSpruce + MeanDaysAboveZero + GroupSize + I(GroupSize^2), data = CurCovarsRepro)
			LinPredsRepro <- ReproData %*% MRepro_Coef
			RetJuv <- exp(LinPredsRepro)
			RetJuv <- data.frame(CurCovarsRepro$GridID, rpois(n = length(RetJuv), lambda = RetJuv))
			colnames(RetJuv) <- c("GridID", "X")
# add to new SJLog
			SJLogAut$X[match(RetJuv$GridID, SJLogAut$GridID)] <- RetJuv$X
# add to Tally
			SJTally$RetJuv[y * 4 - 2] <- sum(SJLogAut$X)
		
########## 2) Set up Mark data: Add GSTerr & GSNbh to Mark covariate data

			CurCovarsE <- merge(CurCovarsE, GSTerrAdd, by = "GridID")			
			CurCovarsE <- CurCovarsE[order(CurCovarsE$GridID), ]
			colnames(CurCovarsE)[which(colnames(CurCovarsE) == "GroupSize")] <- "GSTerr"

			CurCovarsC <- CurCovarsC[CurCovarsC$GridID %in% CurSJ$GridID, ]	
			CurCovarsC <- CurCovarsC[order(CurCovarsC$GridID), ]	

			CurCovarsNB3ToB9 <- merge(CurCovarsNB3ToB9, GSNbhAdd, by = "GridID")	
			CurCovarsNB3ToB9 <- CurCovarsNB3ToB9[order(CurCovarsNB3ToB9$GridID), ]

########## 2) Survival 
	
			DataE <- model.matrix(~ stratumE + NbhAge + GSTerr, data = CurCovarsE)
			DataC <- model.matrix(~ NbhAge, data = CurCovarsC)
			
			LinPreds <- DataE %*% BetasE
			Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
			MMPredsE <- data.frame(GridID = CurCovarsE$GridID, Preds)

			LinPreds <- DataC %*% BetasC
			Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
			MMPredsC <- data.frame(GridID = CurCovarsC$GridID, Preds)

####### for each existing individual simulate using the fitted value
# order by GridID
			MMPredsE <- MMPredsE[order(MMPredsE$GridID), ]
			MMPredsC <- MMPredsC[order(MMPredsC$GridID), ]
			SJLogAut <- SJLogAut[order(SJLogAut$GridID), ]
			CurSJ <- CurSJ[order(CurSJ$GridID), ]

# store the mean predicted survival probabilities (calculated across all individuals to which the model was applied)
			MarkMeanEstimate$S_E[y] <- sum(MMPredsE$Preds[CurSJ$E > 0] * CurSJ$E[CurSJ$E > 0]) / sum(CurSJ$E)
			MarkMeanEstimate$S_C[y] <- sum(MMPredsC$Preds[CurSJ$C > 0] * CurSJ$C[CurSJ$C > 0]) / sum(CurSJ$C)

# draw
			SJLogAut$E <- rbinom(n = nrow(MMPredsE), size = CurSJ$E, prob = MMPredsE$Preds)	
			SJLogAut$C <- rbinom(n = nrow(MMPredsC), size = CurSJ$C, prob = MMPredsC$Preds)
				
			TallyB <- sum(SJLogAut$E)
			TallyNB <- sum(SJLogAut$C)

			SJTally$B[y * 4 - 3] <- sum(SJLogAut$E)
			SJTally$NB[y * 4 - 3] <- sum(SJLogAut$C)	

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterSurv.csv", sep = ""), row.names = F)					
				
########## 3) Transition

# all surviving breeders transition to F

			SJLogAut$F <- SJLogAut$E
			SJLogAut$E <- 0		

# predict using model.matrix and betas

			DataNB3ToB9 <- model.matrix(~ GSNbh + MVol + GSNbh:MVol, data = CurCovarsNB3ToB9)
			LinPreds <- DataNB3ToB9 %*% BetasNB3ToB9
			Preds <- exp(LinPreds) / (1 + exp(LinPreds) + 4)						
			MMPredsNB3ToB9 <- data.frame(GridID = CurCovarsNB3ToB9$GridID, Preds)	

# order by GridID
			MMPredsNB3ToB9 <- MMPredsNB3ToB9[order(MMPredsNB3ToB9$GridID), ]
			SJLogAut <- SJLogAut[order(SJLogAut$GridID), ]
			CurSJ <- CurSJ[order(CurSJ$GridID), ]

# store the mean predicted transition probabilities (calculated across all individuals to which the model was applied)
			MarkMeanEstimate$Psi_CF[y] <- sum(MMPredsNB3ToB9$Preds[CurSJ$C > 0] * CurSJ$C[CurSJ$C > 0]) / sum(CurSJ$C)

# for each existing non-breeder simulate using the fitted value												
			SJLogAut$CToF <- 0
			for (i in 1 : nrow(SJLogAut)){
				if(SJLogAut$C[i] > 0){
					SimMulti <- rmultinom(1, size = SJLogAut$C[i], prob = c(MMPredsNB3ToB9$Preds[i], 1 - MMPredsNB3ToB9$Preds[i]))
					SJLogAut$CToF[i] <- SimMulti[1,]			
					SJLogAut$D[i] <- SimMulti[2,]
				}
			}
			SJLogAut$C <- 0

##### Those that transition C to F will be added to F during Emigration process below (because the emigration model is specific for BToB, NBToB, etc.. For now they are stored in column CToF

			TallyB <- TallyB + sum(SJLogAut$CToF)
			TallyNB <- sum(SJLogAut$D)

			SJTally$TransNB[y * 4 - 2] <- sum(SJLogAut$CToF)

write.csv(SJLogAut, paste(TempFolder,"\\TempAutAfterTrans.csv", sep = ""), row.names = F)		

########## 4) Dispersed juveniles randomly allocated (for their number I follow Kate who used c = 2 x 0.8)

			DispJuv <- round(2 * 0.8 * sum(SJLogAut$X))
			Sel <- sample(1 : nrow(SJLogAut), DispJuv, replace = T)
			Sel <- table(Sel)
			SJLogAut$Z[as.numeric(names(Sel))] <- Sel

			SJTally$DispJuv[y * 4 - 2] <- sum(SJLogAut$Z)

########## 5) Identify existing breeders that emigrate (this is the spring emigration model (breeder = E) as they emigrate between March and Sept

# for breeders
			CurB <- SJLogAut[SJLogAut$F > 0, ]
# the covars and order by GridID
			CurCovarsEmiBToBMarch <- CurCovarsEmiBToBMarch[CurCovarsEmiBToBMarch$GridID %in% CurB$GridID, ]
			CurCovarsEmiBToBMarch <- CurCovarsEmiBToBMarch[order(CurCovarsEmiBToBMarch$GridID), ]
# simulate from model
			EmiData <- model.matrix(~ stg3F + PatchVAMeanAge10km + WinterTemp + WinterPrec + SpringPrec, data = CurCovarsEmiBToBMarch)
			colnames(EmiData)[4:6] <- c("B9:Temp", "B9:Prec", "Prec:B3")
			LinPredsEmi <- EmiData %*% MEmiBToB_Coef
			PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
			PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurB$F, prob = PredsEmi)
			PredsEmi <- data.frame(CurCovarsEmiBToBMarch$GridID, PredsEmi)
			colnames(PredsEmi) <- c("GridID", "EmiBToB")
# add to new SJLog
			SJLogAut$EmiBToB <- 0
			SJLogAut$EmiBToB[match(PredsEmi$GridID, SJLogAut$GridID)] <- PredsEmi$EmiBToB
# and remove from breeder log
			SJLogAut$F <- SJLogAut$F - SJLogAut$EmiBToB

# force breeder to emigrate if there are more than 2 B
			if(nrow(SJLogAut[SJLogAut$F > 2, ]) > 0){
				SJLogAut[SJLogAut$F > 2, ]$EmiBToB <- SJLogAut[SJLogAut$F > 2, ]$EmiBToB + SJLogAut[SJLogAut$F > 2, ]$F - 2	
				SJLogAut[SJLogAut$F > 2, ]$F <- 2
			}

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterEmiOldB.csv", sep = ""), row.names = F)	

########## 6) Identify new breeders (transitioned from non-breeders) that emigrate

			SJLogAut$EmiNBToB <- 0
# if there are new breeders
			CurNBToB <- SJLogAut[SJLogAut$CToF > 0, ]
			if(nrow(CurNBToB) > 0) {
				CurCovarsEmiNBToB <- data.frame(stg3 = factor("C", levels = c("C", "D", "X", "Z")), NoRet = factor(CurNBToB$F, levels = c("0", "1", "2")))		
# simulate from model
				EmiData <- model.matrix(~ stg3 + NoRet, data = CurCovarsEmiNBToB)
				LinPredsEmi <- EmiData %*% MEmiNBToB_Coef
				PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
				PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurNBToB$CToF, prob = PredsEmi)
				PredsEmi <- data.frame(CurNBToB$GridID, PredsEmi)
				colnames(PredsEmi) <- c("GridID", "EmiNBToB")
# add to new SJLog
				SJLogAut$EmiNBToB[match(PredsEmi$GridID, SJLogAut$GridID)] <- PredsEmi$EmiNBToB
			}

# remove those that emigrate from the CToF count
			SJLogAut$CToF <- SJLogAut$CToF - SJLogAut$EmiNBToB
# add the remaining new breeders to F and remove column CToF
			SJLogAut$F <- SJLogAut$F + SJLogAut$CToF
			SJLogAut$CToF <- NULL

# if there are more than 2 B, force the 'new' breeders to emigrate
			if(nrow(SJLogAut[SJLogAut$F > 2, ]) > 0){
				SJLogAut[SJLogAut$F > 2, ]$EmiNBToB <- SJLogAut[SJLogAut$F > 2, ]$EmiNBToB + SJLogAut[SJLogAut$F > 2, ]$F - 2	
				SJLogAut[SJLogAut$F > 2, ]$F <- 2
			}
			
# update SJTally
			SJTally$EmiB[y * 4 - 2] <- sum(SJLogAut$EmiBToB) + sum(SJLogAut$EmiNBToB)

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterEmiNewB.csv", sep = ""), row.names = F)

########## 7) for EmiBToB: if there is a breeder vacancy, let the emigrant immigrate if there is suitable habitat the emigrant can reach, 
########## otherwise move the emigrant to the colonizer column

# predict habitat suitability for each cell
			HSData <- model.matrix(~ PercVAOF + MeanAge + PercVASYF + TempJF + PrecAM + DTM + I(DTM^2) + TempJF:PrecAM, data = CurCovarsHS)
			LinPredsHS <- HSData %*% SF3_Coef
			PredsHS <- exp(LinPredsHS) / (1 + exp(LinPredsHS) )
			PredsHS <- data.frame(CurCovarsHS[, c("GridID", "x", "y")], PredsHS)
			colnames(PredsHS) <- c("GridID", "x", "y",  "HS")
# column number for Sept breeder (F)
			BInd <- 8						
# a column for those that have to colonize
			SJLogAut$ColB <- 0
# a column for cells from which experienced breeders emigrate and which will not be available for immigration by experienced breeders
			SJLogAut$NoVac <- 0
			if(sum(SJLogAut$EmiBToB) > 0) {SJLogAut[SJLogAut$EmiBToB > 0, ]$NoVac <- 1}
			
# immigration at breeder vacancies....
			while(sum(SJLogAut$EmiBToB) > 0){
# identify emigrating birds
				CurEmi <- SJLogAut[SJLogAut$EmiBToB > 0, ]
# let individual with worst habitat suitability disperse first
				CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
				CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
				CurEmi <- CurEmi[1, ]
				SucImm <- FunBImm(CurEmi, SJLogAut, PredsHS, BInd)
# if there is successful immigration
				if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
					SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] <- SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiBToB <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiBToB - 1
				} else {
# ... else move the emigrant to the colonizer column
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$ColB <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$ColB + 1
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiBToB <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiBToB - 1
				}
			}
# remove EmiBToB column
			SJLogAut$EmiBToB <- NULL
		
#################### EmiBToB: for the remaining emigrants, let one colonize a cell with at least the same habitat suitability if possible by distance, let the next emigrant immigrate if possible, 
# otherwise colonize. Those that don't disperse have to stay in same cell

# a column for colonizers that failed to colonize (and will try again to immigrate once all colonizers have gone)
			SJLogAut$FailedColB <- 0

			if(sum(SJLogAut$ColB) > 0){
				CurCols <- SJLogAut[SJLogAut$ColB > 0, ]
				for (i in 1 : sum(CurCols$ColB)){
					CurCol <- SJLogAut[SJLogAut$ColB > 0, ]
# add HS value and order by habitat
					CurCol <- PredsHS[PredsHS$GridID %in% CurCol$GridID, ]
					CurCol <- CurCol[order(CurCol$HS, decreasing = F), ]
# select emigrant in worst habitat
					CurCol <- CurCol[1, ]
# if there was already a colonizer, try immigration first
					if(i > 1){		
						SucImm <- FunBImm(CurCol, SJLogAut, PredsHS, BInd)
						if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
							SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] <- SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
							SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB <- SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB - 1
# if immigration successful, set CurCol to Null
							CurCol <- NULL
						}
					}	# end try immigration
# if this is the first of the emigrants, or if immigration was not successful, try colonization
					if(length(CurCol) > 0){		
						ColSuc <- FunBCol(SJLogAut, PredsHS, BInd, CurCol)
			
# if a 1 was drawn, let the emigrant colonize at the chosen grid cell
						if(length(ColSuc) > 0){
							AddCol <- SJLogAut[1, ]
							AddCol[, c("GridID", "x", "y")] <- c(ColSuc$GridID, ColSuc$x, ColSuc$y)
							Inds <- match("E", colnames(AddCol))
							AddCol[ , Inds : ncol(AddCol)] <- 0
							AddCol[, BInd] <- 1
							SJLogAut <- rbind(SJLogAut, AddCol)
# remove the breeder from the list of colonizers
							SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB <- SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB - 1
# sort by GridID
							SJLogAut <- SJLogAut[order(SJLogAut$GridID), ]
							if(any(duplicated(SJLogAut$GridId))) stop("Duplicated grid numbers")
						} else {
# move the colonizer to the non-colonizing emigrants
							SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB <- SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$ColB - 1
							SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$FailedColB <- SJLogAut[SJLogAut$GridID == CurCol$GridID, ]$FailedColB + 1
						}
					}	# end if try colonization
				}		# end for loop
			}			# end if there are colonizers
# if there are any breeders that couldn't disperse, add them to the breeders in the cell of origin
			SJLogAut$F <- SJLogAut$F + SJLogAut$FailedColB				
# remove redundant columns
			SJLogAut$ColB <- NULL
			SJLogAut$FailedColB <- NULL
			SJLogAut$NoVac <- NULL
# order by GridID
			SJLogAut <- SJLogAut[order(SJLogAut$GridID), ]

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterBImmCol.csv", sep = ""), row.names = F)

########### 8) Let emigrating 'new' breeders immigrate at vacancies, if possible by distance and habitat is within 33 % of suitable habitat. Otherwise, they stay put (as breeder)

# a column for those that have to stay put
			SJLogAut$NoEmi <- 0
			while(sum(SJLogAut$EmiNBToB) > 0){
# identify emigrating birds
				CurEmi <- SJLogAut[SJLogAut$EmiNBToB > 0, ]
# let individual with worst habitat suitability disperse first
				CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
				CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
				CurEmi <- CurEmi[1, ]
				SucImm <- FunNBToBImm(CurEmi, SJLogAut, PredsHS, BInd)
# if there is successful immigration
				if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
					SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] <- SJLogAut[SJLogAut$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiNBToB <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiNBToB - 1
				} else {
# ... else move the emigrant to the failed emigration column
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$NoEmi <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$NoEmi + 1
					SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiNBToB <- SJLogAut[SJLogAut$GridID == CurEmi$GridID, ]$EmiNBToB - 1
				}
			}
# let those that can't emigrate become a breeder in the same cell; there will be forced emigration at the next time step, before breeding begins
			SJLogAut$F <- SJLogAut$F + SJLogAut$NoEmi
			SJLogAut$EmiNBToB <- NULL
			SJLogAut$NoEmi <- NULL

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterNewBImm.csv", sep = ""), row.names = F)

############################################ 9) Identify non-breeders that emigrate between March and September
		
			SJLogAut$EmiNBToNB <- 0
			SJLogAut <- FunEmiNBToNB(SJLogAut)

# add to SJTally
			SJTally$EmiNB[y * 4 - 2] <- sum(SJLogAut$EmiNBToNB)

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterEmiNB.csv", sep = ""), row.names = F)

########## 10) let emigrating non-breeders immigrate to existing groups

			SJLogAut<- FunNBImm(SJLogAut)
# add immigrating NB to D
			SJLogAut$D <- SJLogAut$D + SJLogAut$NewNB
# remove the column in which values were stored
			SJLogAut$EmiNBToNB <- NULL
			SJLogAut$NewNB <- NULL

write.csv(SJLogAut, paste(TempFolder, "\\TempAutAfterImmNB.csv", sep = ""), row.names = F)

# remove unoccupied cells
			SJLogAut <- SJLogAut[apply(SJLogAut[ , c("E", "F", "C", "D", "X", "Z")], 1, sum) > 0, ]			

# update Tally
			SJTally$RetJuv[y * 4 - 2] <- sum(SJLogAut$X)
			SJTally$DispJuv[y * 4 - 2] <- sum(SJLogAut$Z)
			SJTally$B[y * 4 - 2] <- sum(SJLogAut$F)
			SJTally$NB[y * 4 - 2] <- sum(SJLogAut$D)
			SJTally$MeanGSTerr[y * 4 - 2] <- ifelse(nrow(SJLogAut) == 0, 0, mean(FunCalcGSTerr(SJLogAut)$GroupSize))
			SJTally$MaxGSTerr[y * 4 - 2] <- ifelse(nrow(SJLogAut) == 0, 0, max(FunCalcGSTerr(SJLogAut)$GroupSize))
			SJTally$NoTerrs[y * 4 - 2] <- nrow(SJLogAut)
			SJTally$Only1B[y * 4 - 2] <- nrow(SJLogAut[SJLogAut$F ==1, ])

# attach to SJLog
			SJLog <- rbind(SJLog, SJLogAut)
			if(nrow(SJLogAut) == 0) break
			
##############################################################################################################################################################################
############################################################################ September to March ##############################################################################
##############################################################################################################################################################################
			if(max(SJLog$Year) < YearMax){
########## the data in September (SJ and covariates)

				CurSJ <- SJLogAut
				CurSJ <- CurSJ[order(CurSJ$GridID), ]
				
# calculate GSTerr and GSNbh
	
				GSTerrAdd <- FunCalcGSTerr(CurSJ)
				GSNbhAdd <- FunCalcGSNbh(CurSJ)

########## Set up empty log

				SJLogSpr <- CurSJ
				SJLogSpr$Season <- 3
				SJLogSpr$Session <- SJLogSpr$Session + 1
				SJLogSpr$Year <- SJLogSpr$Year + 1
				SJLogSpr[, c("E", "F", "C", "D", "X", "Z")] <- 0

########## 11) Set up Mark data: Add GSTerr & GSNbh to Mark covariate data

				CurCovarsF <- merge(CurCovarsF, GSTerrAdd, by = "GridID")			
				CurCovarsF <- CurCovarsF[order(CurCovarsF$GridID), ]
				colnames(CurCovarsF)[which(colnames(CurCovarsF) == "GroupSize")] <- "GSTerr"

				CurCovarsD <- merge(CurCovarsD, GSTerrAdd, by = "GridID")			
				CurCovarsD <- CurCovarsD[order(CurCovarsD$GridID), ]
				colnames(CurCovarsD)[which(colnames(CurCovarsD) == "GroupSize")] <- "GSTerr"
					
				CurCovarsZ <- merge(CurCovarsZ, GSTerrAdd, by = "GridID")			
				CurCovarsZ <- CurCovarsZ[order(CurCovarsZ$GridID), ]
				colnames(CurCovarsZ)[which(colnames(CurCovarsZ) == "GroupSize")] <- "GSTerr"			

				CurCovarsX <- CurCovarsX[CurCovarsX$GridID %in% CurSJ$GridID, ]		
				CurCovarsX <- CurCovarsX[order(CurCovarsX$GridID), ]	

				CurCovarsNB9ToB3 <- merge(CurCovarsNB9ToB3, GSNbhAdd, by = "GridID")	
				CurCovarsNB9ToB3 <- CurCovarsNB9ToB3[order(CurCovarsNB9ToB3$GridID), ]

				CurCovarsRJToB3 <- merge(CurCovarsRJToB3, GSNbhAdd, by = "GridID")	
				CurCovarsRJToB3 <- CurCovarsRJToB3[order(CurCovarsRJToB3$GridID), ]
			
				CurCovarsDJToB3 <- merge(CurCovarsDJToB3, GSNbhAdd, by = "GridID")	
				CurCovarsDJToB3 <- CurCovarsDJToB3[order(CurCovarsDJToB3$GridID), ]

########## 12) Survival via model matrix

				DataF <- model.matrix(~ stratumF + NbhAge + PD + GSTerr + GSTerr:NbhAge, data = CurCovarsF)
				DataD <- model.matrix(~ stratumD + NbhAge + PD + GSTerr + GSTerr:NbhAge, data = CurCovarsD)
				DataX <- model.matrix(~ stratumX + NbhAge, data = CurCovarsX)
				DataZ <- model.matrix(~ stratumZ + NbhAge + PD + WP + GSTerr + GSTerr:NbhAge + PD:GSTerr, data = CurCovarsZ)

				LinPreds <- DataF %*% BetasF
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
				MMPredsF <- data.frame(GridID = CurCovarsF$GridID, Preds)

				LinPreds <- DataD %*% BetasD
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
				MMPredsD <- data.frame(GridID = CurCovarsD$GridID, Preds)

				LinPreds <- DataX %*% BetasX
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
				MMPredsX <- data.frame(GridID = CurCovarsX$GridID, Preds)

				LinPreds <- DataZ %*% BetasZ
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) )
				MMPredsZ <- data.frame(GridID = CurCovarsZ$GridID, Preds)

####### for each existing individual simulate using the fitted value
# order by GridID
				MMPredsF <- MMPredsF[order(MMPredsF$GridID), ]
				MMPredsD <- MMPredsD[order(MMPredsD$GridID), ]
				MMPredsX <- MMPredsX[order(MMPredsX$GridID), ]
				MMPredsZ <- MMPredsZ[order(MMPredsZ$GridID), ]

				SJLogSpr <- SJLogSpr[order(SJLogSpr$GridID), ]
				CurSJ <- CurSJ[order(CurSJ$GridID), ]

# store the mean predicted survival probabilities (calculated across all individuals to which the model was applied)
				MarkMeanEstimate$S_F[y] <- sum(MMPredsF$Preds[CurSJ$F > 0] * CurSJ$F[CurSJ$F > 0]) / sum(CurSJ$F)
				MarkMeanEstimate$S_D[y] <- sum(MMPredsD$Preds[CurSJ$D > 0] * CurSJ$D[CurSJ$D > 0]) / sum(CurSJ$D)
				MarkMeanEstimate$S_X[y] <- sum(MMPredsX$Preds[CurSJ$X > 0] * CurSJ$X[CurSJ$X > 0]) / sum(CurSJ$X)
				MarkMeanEstimate$S_Z[y] <- sum(MMPredsZ$Preds[CurSJ$Z > 0] * CurSJ$Z[CurSJ$Z > 0]) / sum(CurSJ$Z)

# draw
				SJLogSpr$F <- rbinom(n = nrow(MMPredsF), size = CurSJ$F, prob = MMPredsF$Preds)	
				SJLogSpr$D <- rbinom(n = nrow(MMPredsD), size = CurSJ$D, prob = MMPredsD$Preds)
				SJLogSpr$X <- rbinom(n = nrow(MMPredsX), size = CurSJ$X, prob = MMPredsX$Preds)	
				SJLogSpr$Z <- rbinom(n = nrow(MMPredsZ), size = CurSJ$Z, prob = MMPredsZ$Preds)

				TallyB <- sum(SJLogSpr$F)
				TallyNB <- sum(SJLogSpr[, c("D", "X", "Z")])	

# add to Tally
				SJTally$RetJuv[y * 4 - 1] <- sum(SJLogSpr$X)
				SJTally$DispJuv[y * 4 - 1] <- sum(SJLogSpr$Z)
				SJTally$B[y * 4 - 1] <- sum(SJLogSpr$F)
				SJTally$NB[y * 4 -1] <- sum(SJLogSpr$D)

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterSurv.csv", sep = ""), row.names = F)
			
########## 13) Transition

# all surviving breeders transition to E

				SJLogSpr$E <- SJLogSpr$F
				SJLogSpr$F <- 0
		
# non-breeder
				DataNB9ToB3 <- model.matrix(~ NB9ToB3 + GSNbh + MVol + PD + GSNbh:MVol + PD:GSNbh, data = CurCovarsNB9ToB3)
				LinPreds <- DataNB9ToB3 %*% BetasNB9ToB3
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) + 4)						
				MMPredsNB9ToB3 <- data.frame(GridID = CurCovarsNB9ToB3$GridID, Preds)	
# order by GridID
				MMPredsNB9ToB3 <- MMPredsNB9ToB3[order(MMPredsNB9ToB3$GridID), ]
				SJLogSpr <- SJLogSpr[order(SJLogSpr$GridID), ]
				CurSJ <- CurSJ[order(CurSJ$GridID), ]

# store the mean predicted transition probabilities (calculated across all individuals to which the model was applied)
				MarkMeanEstimate$Psi_DE[y] <- sum(MMPredsNB9ToB3$Preds[CurSJ$D > 0] * CurSJ$D[CurSJ$D > 0]) / sum(CurSJ$D)

# for each existing case simulate using the fitted value												
				SJLogSpr$DToB <- 0
				for (i in 1 : nrow(SJLogSpr)){
					if(SJLogSpr$D[i] > 0){
						SimMulti <- rmultinom(1, size = SJLogSpr$D[i], prob = c(MMPredsNB9ToB3$Preds[i], 1 - MMPredsNB9ToB3$Preds[i]))
						SJLogSpr$DToB[i] <- SimMulti[1,]			
						SJLogSpr$C[i] <- SimMulti[2,]
					}
				}
				SJLogSpr$D <- 0

# retained juvenile
				DataRJToB3 <- model.matrix(~ RJToB3 + GSNbh + MVol + PD + WP + GSNbh:MVol, data = CurCovarsRJToB3)
				LinPreds <- DataRJToB3 %*% BetasRJToB3
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) + 4)						
				MMPredsRJToB3 <- data.frame(GridID = CurCovarsRJToB3$GridID, Preds)	
# order by GridID
				MMPredsRJToB3 <- MMPredsRJToB3[order(MMPredsRJToB3$GridID), ]

# store the mean predicted transition probabilities (calculated across all individuals to which the model was applied)
				MarkMeanEstimate$Psi_XE[y] <- sum(MMPredsRJToB3$Preds[CurSJ$X > 0] * CurSJ$X[CurSJ$X > 0]) / sum(CurSJ$X)

# for each existing case simulate using the fitted value
				SJLogSpr$XToB <- 0
				SJLogSpr$XToNB <- 0											
				for (i in 1 : nrow(SJLogSpr)){
					if(SJLogSpr$X[i] > 0){
						SimMulti <- rmultinom(1, size = SJLogSpr$X[i], prob = c(MMPredsRJToB3$Preds[i], 1 - MMPredsRJToB3$Preds[i]))
						SJLogSpr$XToB[i] <- SimMulti[1,]				
						SJLogSpr$XToNB[i] <- SimMulti[2,]			
					}
				}
				SJLogSpr$X <- 0

# dispersed juvenile: predict using model.matrix and betas

				DataDJToB3 <- model.matrix(~ DJToB3 + GSNbh + MVol + PD + WP + GSNbh:MVol, data = CurCovarsDJToB3)
				LinPreds <- DataDJToB3 %*% BetasDJToB3
				Preds <- exp(LinPreds) / (1 + exp(LinPreds) + 4)						
				MMPredsDJToB3 <- data.frame(GridID = CurCovarsDJToB3$GridID, Preds)	
# order by GridID
				MMPredsDJToB3 <- MMPredsDJToB3[order(MMPredsDJToB3$GridID), ]

# store the mean predicted transition probabilities (calculated across all individuals to which the model was applied)
				MarkMeanEstimate$Psi_ZE[y] <- sum(MMPredsDJToB3$Preds[CurSJ$Z > 0] * CurSJ$Z[CurSJ$Z > 0]) / sum(CurSJ$Z)

# for each existing case simulate using the fitted value	
				SJLogSpr$ZToB <- 0
				SJLogSpr$ZToNB <- 0											
				for (i in 1 : nrow(SJLogSpr)){
					if(SJLogSpr$Z[i] > 0){
						SimMulti <- rmultinom(1, size = SJLogSpr$Z[i], prob = c(MMPredsDJToB3$Preds[i], 1 - MMPredsDJToB3$Preds[i]))
						SJLogSpr$ZToB[i] <- SimMulti[1,]				
						SJLogSpr$ZToNB[i] <- SimMulti[2,]			
					}
				}
				SJLogSpr$Z <- 0

				TallyB <- TallyB + sum(SJLogSpr$DToB + SJLogSpr$XToB + SJLogSpr$ZToB)
				TallyNB <- sum(SJLogSpr$C + SJLogSpr$XToNB + SJLogSpr$ZToNB)

# add to Tally
				SJTally$TransNB[y * 4] <- sum(SJLogSpr$DToB)
				SJTally$TransRJ[y * 4] <- sum(SJLogSpr$XToB)
				SJTally$TransDJ[y * 4] <- sum(SJLogSpr$ZToB)

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterTrans.csv", sep = ""), row.names = F)
					
########## 14) Identify existing breeders that emigrate (this is the autumn emigration model) 

# if there are breeders
				CurB <- SJLogSpr[SJLogSpr$E > 0, ]
# the covars and order by GridID
				CurCovarsEmiBToBSept <- CurCovarsEmiBToBSept[CurCovarsEmiBToBSept$GridID %in% CurB$GridID, ]
				CurCovarsEmiBToBSept <- CurCovarsEmiBToBSept[order(CurCovarsEmiBToBSept$GridID), ]
# simulate from model
				EmiData <- model.matrix(~ stg3F + PatchVAMeanAge10km + WinterTemp + WinterPrec + SpringPrec, data = CurCovarsEmiBToBSept)
				colnames(EmiData)[4:6] <- c("B9:Temp", "B9:Prec", "Prec:B3")
				LinPredsEmi <- EmiData %*% MEmiBToB_Coef
				PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
				PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurB$E, prob = PredsEmi)
				PredsEmi <- data.frame(CurCovarsEmiBToBSept$GridID, PredsEmi)
				colnames(PredsEmi) <- c("GridID", "EmiB")
# add to new SJLog
				SJLogSpr$EmiBToB <- 0
				SJLogSpr$EmiBToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- PredsEmi$EmiB
# remove from the log of breeders
				SJLogSpr$E <- SJLogSpr$E - SJLogSpr$EmiBToB

# if there are more than 2 B, force breeder to emigrate
				if(nrow(SJLogSpr[SJLogSpr$E > 2, ]) > 0){
					SJLogSpr[SJLogSpr$E > 2, ]$EmiBToB <- SJLogSpr[SJLogSpr$E > 2, ]$EmiBToB + SJLogSpr[SJLogSpr$E > 2, ]$E - 2	
					SJLogSpr[SJLogSpr$E > 2, ]$E <- 2
				}

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterEmiOldB.csv", sep = ""), row.names = F)
	

########## 15) Identify 'new' breeders that emigrate
				SJLogSpr$EmiNBToB <- 0	

########## 15.1) For D to breeder
# if there are DToB
				CurToB <- SJLogSpr[SJLogSpr$DToB > 0, ]
				if(nrow(CurToB) > 0) {
					CurCovarsEmiNBToB <- data.frame(stg3 = factor("D", levels = c("C", "D", "X", "Z")), NoRet = factor(CurToB$E, levels = c("0", "1", "2")))		
# simulate from model
					EmiData <- model.matrix(~ stg3 + NoRet, data = CurCovarsEmiNBToB)
					LinPredsEmi <- EmiData %*% MEmiNBToB_Coef
					PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
					PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurToB$DToB, prob = PredsEmi)
					PredsEmi <- data.frame(CurToB$GridID, PredsEmi)
					colnames(PredsEmi) <- c("GridID", "EmiB")
# add to new SJLog
					SJLogSpr$EmiNBToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$EmiNBToB[match(PredsEmi$GridID, SJLogSpr$GridID)] + PredsEmi$EmiB
# remove from DToB log
					SJLogSpr$DToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$DToB[match(PredsEmi$GridID, SJLogSpr$GridID)] - PredsEmi$EmiB
				}

########## 15.2) For X to breeder
				SJLogSpr$EmiJuvToB <- 0
# if there are XToB
				CurToB <- SJLogSpr[SJLogSpr$XToB > 0, ]
				if(nrow(CurToB) > 0) {
					CurCovarsEmiNBToB <- data.frame(stg3 = factor("X", levels = c("C", "D", "X", "Z")), NoRet = factor(CurToB$E, levels = c("0", "1", "2")))
# simulate from model
					EmiData <- model.matrix(~ stg3 + NoRet, data = CurCovarsEmiNBToB)
					LinPredsEmi <- EmiData %*% MEmiNBToB_Coef
					PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
					PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurToB$XToB, prob = PredsEmi)
					PredsEmi <- data.frame(CurToB$GridID, PredsEmi)
					colnames(PredsEmi) <- c("GridID", "EmiB")
# add to new SJLog
					SJLogSpr$EmiJuvToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$EmiJuvToB[match(PredsEmi$GridID, SJLogSpr$GridID)] + PredsEmi$EmiB
# remove from XToB log
					SJLogSpr$XToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$XToB[match(PredsEmi$GridID, SJLogSpr$GridID)] - PredsEmi$EmiB
				}

########## 15.2) For Z to breeder
# if there are ZToB
				CurToB <- SJLogSpr[SJLogSpr$ZToB > 0, ]
				if(nrow(CurToB) > 0) {
					CurCovarsEmiNBToB <- data.frame(stg3 = factor("Z", levels = c("C", "D", "X", "Z")), NoRet = factor(CurToB$E, levels = c("0", "1", "2")))
# simulate from model
					EmiData <- model.matrix(~ stg3 + NoRet, data = CurCovarsEmiNBToB)
					LinPredsEmi <- EmiData %*% MEmiNBToB_Coef
					PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
					PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurToB$ZToB, prob = PredsEmi)
					PredsEmi <- data.frame(CurToB$GridID, PredsEmi)
					colnames(PredsEmi) <- c("GridID", "EmiB")
# add to new SJLog
					SJLogSpr$EmiJuvToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$EmiJuvToB[match(PredsEmi$GridID, SJLogSpr$GridID)] + PredsEmi$EmiB
# remove from ZToB log
					SJLogSpr$ZToB[match(PredsEmi$GridID, SJLogSpr$GridID)] <- SJLogSpr$ZToB[match(PredsEmi$GridID, SJLogSpr$GridID)] - PredsEmi$EmiB
				}

# add non-emigrating new breeders to E and remove the column that logged them
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$DToB
				SJLogSpr$DToB <- NULL
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$XToB
				SJLogSpr$XToB <- NULL
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$ZToB
				SJLogSpr$ZToB <- NULL

# if there are more than 2 B, force breeder to emigrate
				if(nrow(SJLogSpr[SJLogSpr$E > 2, ]) > 0){
					SJLogSpr[SJLogSpr$E > 2, ]$EmiBToB <- SJLogSpr[SJLogSpr$E > 2, ]$EmiBToB + SJLogSpr[SJLogSpr$E > 2, ]$E - 2	
					SJLogSpr[SJLogSpr$E > 2, ]$E <- 2
				}
# add emigrating breeders to Tally
				SJTally$EmiB[y * 4] <- sum(SJLogSpr$EmiBToB) + sum(SJLogSpr$EmiNBToB) + sum(SJLogSpr$EmiJuvToB)

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterEmiNewB.csv", sep = ""), row.names = F)


########## for EmiBToB: if there is a vacancy, let the emigrant immigrate if there is suitable habitat the emigrant can reach, otherwise move the emigrant to the colonizer column

# predict habitat suitability for each cell
				HSData <- model.matrix(~ PercVAOF + MeanAge + PercVASYF + TempJF + PrecAM + DTM + I(DTM^2) + TempJF:PrecAM, data = CurCovarsHS)
				LinPredsHS <- HSData %*% SF3_Coef
				PredsHS <- exp(LinPredsHS) / (1 + exp(LinPredsHS) )
				PredsHS <- data.frame(CurCovarsHS[, c("GridID", "x", "y")], PredsHS)
				colnames(PredsHS) <- c("GridID", "x", "y",  "HS")
# Column index for breeder
				BInd <- 7
# a column for cells from which experienced breeders emigrate and which will not be available for immigration by experienced breeders
				SJLogSpr$NoVac <- 0
				if(sum(SJLogSpr$EmiBToB) > 0) {SJLogSpr[SJLogSpr$EmiBToB > 0, ]$NoVac <- 1}

# a column for those that have to colonize
				SJLogSpr$ColB <- 0

				while(sum(SJLogSpr$EmiBToB) > 0){
# identify emigrating birds
					CurEmi <- SJLogSpr[SJLogSpr$EmiBToB > 0, ]
# let individual with worst habitat suitability disperse first
					CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
					CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
					CurEmi <- CurEmi[1, ]
					SucImm <- FunBImm(CurEmi, SJLogSpr, PredsHS, BInd)
# if there is successful immigration
					if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
						SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiBToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiBToB - 1
					} else {
# ... else move the emigrant to the colonizer column
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$ColB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$ColB + 1
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiBToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiBToB - 1
					}
				}
# remove EmiBToB column
				SJLogSpr$EmiBToB <- NULL
		
#################### EmiBToB: for the remaining emigrants, let one colonize a cell with at least the same habitat suitability if possible by distance, let the next emigrant immigrate if possible, 
# otherwise colonize. Those that don't disperse have to stay in same cell

# a column for colonizers that failed to colonize (and will try again to immigrate once all colonizers have gone)
				SJLogSpr$FailedColB <- 0

				if(sum(SJLogSpr$ColB) > 0){
					CurCols <- SJLogSpr[SJLogSpr$ColB > 0, ]
					for (i in 1 : sum(CurCols$ColB)){
						CurCol <- SJLogSpr[SJLogSpr$ColB > 0, ]
# add HS value and order by habitat
						CurCol <- PredsHS[PredsHS$GridID %in% CurCol$GridID, ]
						CurCol <- CurCol[order(CurCol$HS, decreasing = F), ]
# select emigrant in worst habitat
						CurCol <- CurCol[1, ]
# if there was already a colonizer, try immigration first
						if(i > 1){		
							SucImm <- FunBImm(CurCol, SJLogSpr, PredsHS, BInd)
							if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
								SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
								SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB <- SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB - 1
# if immigration successful, set CurCol to Null
								CurCol <- NULL
							}
						}	# end if try immigration
# if this is the first of the emigrants, or if immigration was not successful, try colonization
						if(length(CurCol) > 0){		
							ColSuc <- FunBCol(SJLogSpr, PredsHS, BInd, CurCol)
			
# if a 1 was drawn, let the emigrant colonize at the chosen grid cell
							if(length(ColSuc) > 0){
								AddCol <- SJLogSpr[1, ]
								AddCol[, c("GridID", "x", "y")] <- c(ColSuc$GridID, ColSuc$x, ColSuc$y)
								Inds <- match("E", colnames(AddCol))
								AddCol[ , Inds : ncol(AddCol)] <- 0
								AddCol[, BInd] <- 1
								SJLogSpr <- rbind(SJLogSpr, AddCol)
# remove the breeder from the list of colonizers
								SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB <- SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB - 1
# sort by GridID
								SJLogSpr <- SJLogSpr[order(SJLogSpr$GridID), ]
								if(any(duplicated(SJLogSpr$GridId))) stop("Duplicated grid numbers")
							} else {
# move the colonizer to the non-colonizing emigrants
								SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB <- SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$ColB - 1
								SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$FailedColB <- SJLogSpr[SJLogSpr$GridID == CurCol$GridID, ]$FailedColB + 1
							}
						}	# end if try colonization
					}		# end for loop
				}			# end if there are colonizers
# if there are any breeders that couldn't disperse, add them to the breeders in the cell of origin
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$FailedColB				
# remove redundant columns
				SJLogSpr$ColB <- NULL
				SJLogSpr$FailedColB <- NULL
				SJLogSpr$NoVac <- NULL
# order by GridID
				SJLogSpr <- SJLogSpr[order(SJLogSpr$GridID), ]

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterBImmCol.csv", sep = ""), row.names = F)

########### Let emigrating new breeders immigrate at vacancies, if possible by distance and habitat is within threshold of suitable habitat. Otherwise, they stay put

# a column for those that have to stay put
				SJLogSpr$NoEmi <- 0
				while(sum(SJLogSpr$EmiNBToB) > 0){
# identify emigrating birds
					CurEmi <- SJLogSpr[SJLogSpr$EmiNBToB > 0, ]
# let individual with worst habitat suitability disperse first
					CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
					CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
					CurEmi <- CurEmi[1, ]
					SucImm <- FunNBToBImm(CurEmi, SJLogSpr, PredsHS, BInd)
# if there is successful immigration
					if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
						SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiNBToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiNBToB - 1
					} else {
# ... else move the emigrant to the non-immigrant column
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoEmi <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoEmi + 1
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiNBToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiNBToB - 1
					}
				}
# let those that can't emigrate become a breeder in the same cell
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$NoEmi
				SJLogSpr$EmiNBToB <- NULL
				SJLogSpr$NoEmi <- NULL

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterNBToBImm.csv", sep = ""), row.names = F)

################ let emmigrating JuvToB immigrate independent of habitat suitability, but dependent on distance

				SJLogSpr$NoEmi <- 0
				while(sum(SJLogSpr$EmiJuvToB) > 0){
# identify emigrating birds
					CurEmi <- SJLogSpr[SJLogSpr$EmiJuvToB > 0, ]
# let individual with worst habitat suitability disperse first
					CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
					CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
					CurEmi <- CurEmi[1, ]
					SucImm <- FunJuvToBImm(CurEmi, SJLogSpr, BInd)
# if there is successful immigration
					if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
						SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the list of emigrants
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiJuvToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiJuvToB - 1
					} else {
# ... else move the emigrant to the no emigration column
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoEmi <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoEmi + 1
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiJuvToB <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$EmiJuvToB - 1
					}
				}
# let those that can't emigrate become a breeder in the same cell
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$NoEmi
				SJLogSpr$EmiJuvToB <- NULL
				SJLogSpr$NoEmi <- NULL

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterJuvToBImm.csv", sep = ""), row.names = F)

################ if there are more than one territory with only one breeder, within the max dispersal distance, force emigration from cell with lowest HS,
################ using dispersal distance probability, until only max one such territory is left
				SJLogSpr$NoVac <- 0
				while(nrow(SJLogSpr[SJLogSpr$E == 1,]) > 1){
					CurEmi <- SJLogSpr[SJLogSpr$E == 1,]
					CurEmi <- PredsHS[PredsHS$GridID %in% CurEmi$GridID, ]
					CurEmi <- CurEmi[order(CurEmi$HS, decreasing = F), ]
					CurEmi <- CurEmi[1, ]
					SucImm <- FunSingleBImm(CurEmi, SJLogSpr, PredsHS, BInd)
# if there is successful immigration
					if(length(SucImm) > 0){ 
# let the emigrating breeder immigrate at the immigration GridID
						SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == SucImm$GridID, BInd] + 1
# remove the breeder from the source cell
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, BInd] - 1
					} else {
# ... else move the emigrant to the failed move column
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoVac <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, ]$NoVac + 1
						SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, BInd] <- SJLogSpr[SJLogSpr$GridID == CurEmi$GridID, BInd] - 1
					}
				}
# if there is no vacancy within the maximum dispersal distance, they have to stay put
				SJLogSpr$E <- SJLogSpr$E + SJLogSpr$NoVac
# remove column that held the breeders that have not moved
				SJLogSpr$NoVac <- NULL
# order by GridID
				SJLogSpr <- SJLogSpr[order(SJLogSpr$GridID), ]

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterSingleBImm.csv", sep = ""), row.names = F)

################ identify emigrating non-breeders
				SJLogSpr$EmiNBToNB <- 0
				SJLogSpr <- FunEmiNBToNB(SJLogSpr)

# add non-emigrating NB to NB tally
				SJLogSpr$C <- SJLogSpr$C + SJLogSpr$XToNB
				SJLogSpr$C <- SJLogSpr$C + SJLogSpr$ZToNB
# remove ZToNB and XToNB columns
				SJLogSpr$XToNB <- NULL
				SJLogSpr$ZToNB <- NULL
# add emigrants to Tally	
				SJTally$EmiNB[y * 4] <- sum(SJLogSpr$EmiNBToNB)

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterBImmi.csv", sep = ""), row.names = F)

########## 10) let emigrating non-breeders immigrate to existing groups

				SJLogSpr<- FunNBImm(SJLogSpr)

# add immigrating NB to C
				SJLogSpr$C <- SJLogSpr$C + SJLogSpr$NewNB
# remove the column in which values were stored
				SJLogSpr$EmiNBToNB <- NULL
				SJLogSpr$NewNB <- NULL

write.csv(SJLogSpr, paste(TempFolder, "\\TempSprAfterNBImmi.csv", sep = ""), row.names = F)

# remove unoccupied cells
				SJLogSpr <- SJLogSpr[apply(SJLogSpr[ , c("E", "F", "C", "D", "X", "Z")], 1, sum) > 0, ]
# add totals to Tally
				SJTally$B[y * 4] <- sum(SJLogSpr$E)
				SJTally$NB[y * 4] <- sum(SJLogSpr$C) 
				SJTally$MeanGSTerr[y * 4] <- ifelse(nrow(SJLogSpr) == 0, 0, mean(FunCalcGSTerr(SJLogSpr)$GroupSize))
				SJTally$MaxGSTerr[y * 4] <- ifelse(nrow(SJLogSpr) == 0, 0, max(FunCalcGSTerr(SJLogSpr)$GroupSize))
				SJTally$NoTerrs[y * 4] <- nrow(SJLogSpr)
				SJTally$Only1B[y * 4] <- nrow(SJLogSpr[SJLogSpr$E == 1, ])			
# attach to SJLog
				SJLog <- rbind(SJLog, SJLogSpr)
			}		# end if											
		}			# end year
		write.csv(SJLog, paste(OutFolder, "\\", paste(TargetAreaName, "SJLog", NamingComment, Rep, ".csv", sep = "_"), sep = ""), row.names = F)
		write.csv(MarkMeanEstimate, paste(OutFolder, "\\", paste(TargetAreaName, "MarkEstimates", NamingComment, Rep, ".csv", sep = "_"), sep = ""), row.names = F)
		SJTally <- rbind(SJTallyInit, SJTally)
		write.csv(SJTally, paste(OutFolder, "\\", paste(TargetAreaName, "Tally", NamingComment, Rep, ".csv", sep = "_"), sep = ""), row.names = F)
	}				# end Rep
}


################################################## Functions to calculate group size at the scale of the territory and neigbhourhood ###############################
####################################################################################################################################################################

# Function to calculate GSTerr on an SJ log (e.g. CurSJ / SJLogAut)

FunCalcGSTerr <- function(MyLog){
	GSTerrAdd <- data.frame(MyLog$GridID, apply(MyLog[, c("E", "F", "C", "D", "X", "Z")], 1, sum))
	colnames(GSTerrAdd) <- c("GridID", "GroupSize")
	return(GSTerrAdd)
}

# Function to calcualte GSNbh on an SJ log (e.g. CurSJ / SJLogAut)

FunCalcGSNbh <- function(MyLog){
	MyCoords <- cbind(MyLog$x, MyLog$y)
	CurSJSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(MyLog), proj4string = CRS("+init=epsg:3021"))

	GSNbhAdd <- as.data.frame(matrix(ncol = 2, nrow = nrow(MyLog)))
	colnames(GSNbhAdd) <- c("GridID", "GSNbh")
	GSNbhAdd[, 1] <- MyLog$GridID
	for (p in 1 : nrow(MyLog)) {
		Obs <- MyLog[p, ]
		Cell10km <- as(extent(Obs$x - 5000, Obs$x + 5000, Obs$y - 5000, Obs$y + 5000), "SpatialPolygons")
		proj4string(Cell10km) <- CRS("+init=epsg:3021")
		SJIn10km <- CurSJSPDF[Cell10km, ]
		SJIn10km <- data.frame(SJIn10km)
		GSNbhAdd[p, 2] <- mean(apply(SJIn10km[, c("E", "F", "C", "D", "X", "Z")], 1, sum))	
	}
	return(GSNbhAdd)
}


############################################################# Functions for experienced breeder immigration at vacancies ###########################################
####################################################################################################################################################################

FunBImm <- function(CurEmi, SJL, PredsHS, BInd){

	SucImm <- NULL

# make SPDF of location of emigrating bird for the single breeder in the cell with lowest HS
	MyCoords <- cbind(CurEmi$x, CurEmi$y)
	CurEmiSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurEmi), proj4string = CRS("+init=epsg:3021"))			
# make SPDF of locations with breeder vacancy...
	CurVacant <- SJL[SJL[BInd] < 2 & SJL$NoVac == 0, ]
	CurVacant <- PredsHS[PredsHS$GridID %in% CurVacant$GridID, ]
# ... if habitat suitability is at least as high as in the emigrating location...
	CurVacant <- CurVacant[CurVacant$HS >= CurEmi$HS, ]
	if(nrow(CurVacant) > 0){
		MyCoords <- cbind(CurVacant$x, CurVacant$y)
		CurVacantSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurVacant), proj4string = CRS("+init=epsg:3021"))
# calculate distances between emigrant and vacancies
		Dists <- gDistance(CurEmiSPDF, CurVacantSPDF, byid = T)
		Dists <- data.frame(GridID = CurVacant$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance and closer than min observed dispersal distance	
		Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# add the fitted frequencies from the dispersal distance bin model
		if(nrow(Dists) > 0) {
			DistData <- model.matrix(~ Dists$Dist)
			LinPredsDist <- DistData %*% DispDistCoef
			Dists$FittedFreq <- exp(LinPredsDist)
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
# observed frequency of dispersal events in first 1 km bin was 128.
			if(min(Dists$Dist) <=500){
				Dists[Dists$Dist <= 500, ]$FittedFreq <- (Dists[Dists$Dist <= 500, ]$Dist - 278) * 128/222
			}
# convert frequency to probability
			Dists$DispProb <- Dists$FittedFreq / 230			# 230 observed dispersal events
# if there is a close vacancy, immigration will always take place, otherwise only if a draw from the dispersal probabilities results in at least one 1
# as close I use (5000000/pi)^0.5, because the upper limit of reported home range sizes is 5 sqkm. I assume that within this area Sibe jays will notice if a vacancy occurs.
# This distance (1261.566) is also very similar to the dispersal distance in Griesser et al. 2014 for delayed dispersers (1250 m).
			DispSuc <- NULL
			if(min(Dists$Dist) < (5000000/pi)^0.5){DispSuc <- 1} else {DispSuc <- rbinom(nrow(Dists), 1, prob = Dists$DispProb)}
# if at least one of the draws results in a 1, let emigrant immigrate using Prob(HS) * Prob(Disp)
			if(sum(DispSuc) > 0){
				Dists <- merge(Dists, PredsHS, by = "GridID")
				Dists$CombProb <- Dists$DispProb * Dists$HS
				Dists <- Dists[order(Dists$CombProb, decreasing = T), ]
				Dists$CumSum <- cumsum(Dists$CombProb)	
				Draw <- runif(1, min = 0, max = sum(Dists$CombProb))
				SucImm <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
			}
		}
	}
	return(SucImm)
}

############################################################### Functions for breeder colonization at unoccupied cells #############################################
####################################################################################################################################################################

FunBCol <- function(SJL, PredsHS, BInd, CurCol){
	ColSuc <- NULL	
	MyCoords <- cbind(CurCol$x, CurCol$y)
	CurColSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurCol), proj4string = CRS("+init=epsg:3021"))
# identify empty cells
	EmptyCells <- PredsHS[!PredsHS$GridID %in% SJL$GridID, ]
# remove cells where the habitat is worse
	EmptyCells <- EmptyCells[EmptyCells$HS >= CurCol$HS, ]		
# if there are cells ...
	if(nrow(EmptyCells) > 0){
# ...make SPDF of empty cells
		MyCoords <- cbind(EmptyCells$x, EmptyCells$y)
		EmptyCellsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(EmptyCells), proj4string = CRS("+init=epsg:3021"))
# calculate distances between emigrant and empty cells
		Dists <- gDistance(CurColSPDF, EmptyCellsSPDF, byid = T)
		Dists <- data.frame(GridID = EmptyCellsSPDF$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance	
		Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# if there are any empty cells within max observed dispersal distance...
		if(nrow(Dists) > 0){
# remove vacancies that are too close to an already occupied territory
# ...make SPDF of remaining empty cells
			EmptyCells <- EmptyCells[EmptyCells$GridID %in% Dists$GridID, ]
			MyCoords <- cbind(EmptyCells$x, EmptyCells$y)
			EmptyCellsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(EmptyCells), proj4string = CRS("+init=epsg:3021"))
# ...make SPDF of occupied cells
			MyCoords <- cbind(SJL$x, SJL$y)
			OccCellsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(SJL), proj4string = CRS("+init=epsg:3021"))
# calculate distances between empty and occupied cells
			Dists <- gDistance(EmptyCellsSPDF, OccCellsSPDF, byid = T)
			MinDists <- numeric(length = nrow(EmptyCells))
			for (i in 1 : ncol(Dists)){
				MinDists[i] <- min(Dists[, i])
			}
			EmptyToOccDists <- data.frame(GridID = EmptyCellsSPDF$GridID, Dist = MinDists)
# calculate distances between occupied cells
			Dists <- gDistance(OccCellsSPDF, byid = T)
			MinDists <- numeric(length = nrow(Dists))
			for (i in 1 : nrow(Dists)){
				CurDists <- Dists[i, ]
				CurDists <- CurDists[order(CurDists)]
				MinDists[i] <- CurDists[2]
			}
			BetweenOccDists <- data.frame(GridID = OccCellsSPDF$GridID, Dist = MinDists)
# calculate how many should be in the distance bins up to 750 m based on the observed data
			brk <- seq(-5, 50245, by = 250)			# this puts 250 m in the Max500 distance bin
			Hist <- hist(BetweenOccDists$Dist, breaks = brk)
			TargetNoMax500 <- round(nrow(SJL) * ObsTerrDistProp[2])
			TargetNoMax750 <- round(nrow(SJL) * ObsTerrDistProp[3])
# are any more allowed in these distance bins?
			TargetNoMax500 <- TargetNoMax500 - Hist$counts[2]
			TargetNoMax750 <- TargetNoMax750 - Hist$counts[3]
# if none are allowed, remove empty cells in these distance bins
			if(TargetNoMax500 <= 0){EmptyToOccDists <- EmptyToOccDists[EmptyToOccDists$Dist >= 500, ]}
			if(TargetNoMax750 <= 0){EmptyToOccDists <- EmptyToOccDists[EmptyToOccDists$Dist < 500 | EmptyToOccDists$Dist >= 750, ]}
# add the fitted frequencies from the dispersal distance bin model
			if(nrow(EmptyToOccDists) > 0){
				DistData <- model.matrix(~ EmptyToOccDists$Dist)
				LinPredsDist <- DistData %*% DispDistCoef
				EmptyToOccDists$FittedFreq <- as.vector(exp(LinPredsDist))
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
# observed frequency of dispersal events in first 1 km bin was 128.
				if(min(EmptyToOccDists$Dist) <=500){
					EmptyToOccDists[EmptyToOccDists$Dist <= 500, ]$FittedFreq <- (EmptyToOccDists[EmptyToOccDists$Dist <= 500, ]$Dist - 278) * 128/222
				}
				EmptyToOccDists$DispProb <- EmptyToOccDists$FittedFreq / 230	# convert the frequency with which dispersal occurs at this distance to probability (a total of 230 dispersal events were observed)
# change negative to 0 probability
				if(min(EmptyToOccDists$DispProb) < 0){EmptyToOccDists[EmptyToOccDists$DispProb < 0, ]$DispProb <- 0}
# draw from dispersal probability
				DispSuc <- NULL
				DispSuc <- rbinom(nrow(EmptyToOccDists), 1, prob = EmptyToOccDists$DispProb)
# if at least one the draws results in a 1, let emigrant immigrate using Prob(HS) * Prob(Disp)
				if(sum(DispSuc) > 0){
					Dists <- merge(EmptyToOccDists, PredsHS, by = "GridID")
					Dists$CombProb <- Dists$DispProb * Dists$HS
					Dists <- Dists[order(Dists$CombProb, decreasing = T), ]
					Dists$CumSum <- cumsum(Dists$CombProb)						
					Draw <- runif(1, min = 0, max = sum(Dists$CombProb))
					ColSuc <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
				}
			}
		}
	}
	return(ColSuc)
}	

################################################################ Functions for 'new' breeder immigration at vacancies ##############################################
####################################################################################################################################################################

FunNBToBImm <- function(CurEmi, SJL, PredsHS, BInd){

	SucImm <- NULL

# make SPDF of location of emigrating bird for the single breeder in the cell with lowest HS
	MyCoords <- cbind(CurEmi$x, CurEmi$y)
	CurEmiSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurEmi), proj4string = CRS("+init=epsg:3021"))			
# make SPDF of locations with breeder vacancy...
	CurVacant <- SJL[SJL$GridID != CurEmi$GridID & SJL[BInd] < 2, ]
	CurVacant <- PredsHS[PredsHS$GridID %in% CurVacant$GridID, ]
# habitat suitability within 33% of origin
	ThreshHS <- CurEmi$HS - (CurEmi$HS * 0.33)
# ... if habitat suitability is at least as high as the threshold habitat suitability...
	CurVacant <- CurVacant[CurVacant$HS >= ThreshHS, ]
	if(nrow(CurVacant) > 0){
		MyCoords <- cbind(CurVacant$x, CurVacant$y)
		CurVacantSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurVacant), proj4string = CRS("+init=epsg:3021"))
# calculate distances between emigrant and vacancies
		Dists <- gDistance(CurEmiSPDF, CurVacantSPDF, byid = T)
		Dists <- data.frame(GridID = CurVacant$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance	
		Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# add the fitted frequencies from the dispersal distance bin model
		if(nrow(Dists) > 0){
			DistData <- model.matrix(~ Dists$Dist)
			LinPredsDist <- DistData %*% DispDistCoef
			Dists$FittedFreq <- exp(LinPredsDist)
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
# observed frequency of dispersal events in first 1 km bin was 128.
			if(min(Dists$Dist) <=500){
				Dists[Dists$Dist <= 500, ]$FittedFreq <- (Dists[Dists$Dist <= 500, ]$Dist - 278) * 128/222
			}
# convert frequency to probability
			Dists$DispProb <- Dists$FittedFreq / 230			# 230 observed dispersal events
# if there is a close vacancy, immigration will always take place, otherwise only if a draw from the dispersal probabilities results in at least one 1
# as close I use (5000000/pi)^0.5, because the upper limit of reported home range sizes is 5 sqkm. I assume that within this area Sibe jays will notice if a vacancy occurs.
# This distance (1261.566) is also very similar to the dispersal distance in Griesser et al. 2014 for delayed dispersers (1250 m).
			DispSuc <- NULL
			if(min(Dists$Dist) < (5000000/pi)^0.5){DispSuc <- 1} else {DispSuc <- rbinom(nrow(Dists), 1, prob = Dists$DispProb)}
# if at least one of the draws results in a 1, let emigrant immigrate using Prob(HS) * Prob(Disp)
			if(sum(DispSuc) > 0){
				Dists <- merge(Dists, PredsHS, by = "GridID")
				Dists$CombProb <- Dists$DispProb * Dists$HS
				Dists <- Dists[order(Dists$CombProb, decreasing = T), ]
				Dists$CumSum <- cumsum(Dists$CombProb)	
				Draw <- runif(1, min = 0, max = sum(Dists$CombProb))
				SucImm <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
			}
		}
	}
	return(SucImm)
}

############################################################# Functions for juvenile to breeder immigration at vacancies ###########################################
####################################################################################################################################################################

FunJuvToBImm <- function(CurEmi, SJL, BInd){

	SucImm <- NULL

# make SPDF of location of emigrating bird for the single breeder in the cell with lowest HS
	MyCoords <- cbind(CurEmi$x, CurEmi$y)
	CurEmiSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurEmi), proj4string = CRS("+init=epsg:3021"))			
# make SPDF of locations with breeder vacancy...
	CurVacant <- SJL[SJL$GridID != CurEmi$GridID & SJL[BInd] < 2, ]
	if(nrow(CurVacant) > 0){
		MyCoords <- cbind(CurVacant$x, CurVacant$y)
		CurVacantSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurVacant), proj4string = CRS("+init=epsg:3021"))
# calculate distances between emigrant and vacancies
		Dists <- gDistance(CurEmiSPDF, CurVacantSPDF, byid = T)
		Dists <- data.frame(GridID = CurVacant$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance	
		Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# add the fitted frequencies from the dispersal distance bin model
		if(nrow(Dists) > 0){
			DistData <- model.matrix(~ Dists$Dist)
			LinPredsDist <- DistData %*% DispDistCoef
			Dists$FittedFreq <- exp(LinPredsDist)
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
# observed frequency of dispersal events in first 1 km bin was 128.
			if(min(Dists$Dist) <=500){
				Dists[Dists$Dist <= 500, ]$FittedFreq <- (Dists[Dists$Dist <= 500, ]$Dist - 278) * 128/222
			}
# convert frequency to probability
			Dists$DispProb <- Dists$FittedFreq / 230			# 230 observed dispersal events
# if there is a close vacancy, immigration will always take place, otherwise only if a draw from the dispersal probabilities results in at least one 1
# as close I use (5000000/pi)^0.5, because the upper limit of reported home range sizes is 5 sqkm. I assume that within this area Sibe jays will notice if a vacancy occurs.
# This distance (1261.566) is also very similar to the dispersal distance in Griesser et al. 2014 for delayed dispersers (1250 m).
			DispSuc <- NULL
			if(min(Dists$Dist) < (5000000/pi)^0.5){DispSuc <- 1} else {DispSuc <- rbinom(nrow(Dists), 1, prob = Dists$DispProb)}
# if at least one of the draws results in a 1, let emigrant immigrate using Prob(Disp)
			if(sum(DispSuc) > 0){
				Dists <- Dists[order(Dists$DispProb, decreasing = T), ]
				Dists$CumSum <- cumsum(Dists$DispProb)	
				Draw <- runif(1, min = 0, max = sum(Dists$DispProb))
				SucImm <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
			}
		}
	}
	return(SucImm)
}

################################################################ Function for forced dispersal of single breeders ##################################################
####################################################################################################################################################################

FunSingleBImm <- function(CurEmi, SJL, PredsHS, BInd){

	SucImm <- NULL

# make SPDF of location of emigrating bird for the single breeder in the cell with lowest HS
	MyCoords <- cbind(CurEmi$x, CurEmi$y)
	CurEmiSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurEmi), proj4string = CRS("+init=epsg:3021"))			
# make SPDF of locations with breeder vacancy...
	CurVacant <- SJL[SJL$GridID != CurEmi$GridID & SJL[BInd]== 1, ]
	CurVacant <- PredsHS[PredsHS$GridID %in% CurVacant$GridID, ]
	MyCoords <- cbind(CurVacant$x, CurVacant$y)
	CurVacantSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurVacant), proj4string = CRS("+init=epsg:3021"))
# calculate distances between emigrant and vacancies
	Dists <- gDistance(CurEmiSPDF, CurVacantSPDF, byid = T)
	Dists <- data.frame(GridID = CurVacant$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance	
	Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# add the fitted frequencies from the dispersal distance bin model
	if(nrow(Dists) > 0){
		DistData <- model.matrix(~ Dists$Dist)
		LinPredsDist <- DistData %*% DispDistCoef
		Dists$FittedFreq <- exp(LinPredsDist)
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
# observed frequency of dispersal events in first 1 km bin was 128. 
		if(min(Dists$Dist) <=500){
			Dists[Dists$Dist <= 500, ]$FittedFreq <- (Dists[Dists$Dist <= 500, ]$Dist - 278) * 128/222
		}
# convert frequency to probability
		Dists$DispProb <- Dists$FittedFreq / 230			# 230 observed dispersal events
# single breeders will always disperse
		DispSuc <- NULL
		Dists <- merge(Dists, PredsHS, by = "GridID")
		Dists$CombProb <- Dists$DispProb * Dists$HS
		Dists <- Dists[order(Dists$CombProb, decreasing = T), ]
		Dists$CumSum <- cumsum(Dists$CombProb)	
		Draw <- runif(1, min = 0, max = sum(Dists$CombProb))
		SucImm <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
	}
	return(SucImm)
}

###################################################################### Function for immigration of non-breeders ####################################################
####################################################################################################################################################################

FunNBImm <- function(SJL){

	DrawSuc <- NULL

# add column NewNB to hold the NB's that immigrate in this function
	SJL$NewNB <- 0
# make SPDF of locations of existing groups, with >= 1 breeder
	CurGroups <- SJL[apply(SJL[ , c("E", "F")], 1, sum) > 0, ]
# if there are no existing groups with breeders, the immigrants stays put
	if(nrow(CurGroups) == 0){
		SJL$NewNB <- SJL$NewNB + SJL$EmiNBToNB
		SJL$EmiNBToNB <- 0
	}else{
# the emigrating birds
		NBEmiGroups <- SJL[SJL$EmiNBToNB > 0, ]
		if(nrow(CurGroups) > 0){
			MyCoords <- cbind(CurGroups$x, CurGroups$y)				
			CurGroupsSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurGroups), proj4string = CRS("+init=epsg:3021"))
			while(sum(SJL$EmiNBToNB) > 0){
# make SPDF of location of emigrating bird
				CurEmi <- SJL[SJL$EmiNBToNB > 0, ]
				CurEmi <- CurEmi[1, ]
				MyCoords <- cbind(CurEmi$x, CurEmi$y)
				CurEmiSPDF <- SpatialPointsDataFrame(MyCoords, data = data.frame(CurEmi), proj4string = CRS("+init=epsg:3021"))
# calculate distances
				Dists <- gDistance(CurEmiSPDF, CurGroupsSPDF, byid = T)
				Dists <- data.frame(GridID = CurGroups$GridID, Dist = Dists[, 1])
# remove vacancies further away than max observed dispersal distance	
				Dists <- Dists[Dists$Dist < 13321 & Dists$Dist > 278, ]
# add the fitted frequencies from the dispersal distance bin model
				if(nrow(Dists) > 0){
					DistData <- model.matrix(~ Dists$Dist)
					LinPredsDist <- DistData %*% DispDistCoef
					Dists$FittedFreq <- exp(LinPredsDist)
# for distances within 500, fit a simple line (with slope 128/ 222 (222: 500 - 278 (min obs disp dist))
					if(min(Dists$Dist) <=500){
						Dists[Dists$Dist <= 500, ]$FittedFreq <- (Dists[Dists$Dist <= 500, ]$Dist - 278) * 128/222
					}
					Dists$DispProb <- Dists$FittedFreq / 230	# convert the frequency with which dispersal occurs at this distance to probability (a total of 230 dispersal events were observed)			
# immigration will always take place if a draw from the dispersal probabilities results in at least one 1
					DispSuc <- NULL
					DispSuc <- rbinom(nrow(Dists), 1, prob = Dists$DispProb)
# if at least one of the draws results in a 1, let emigrant immigrate using Prob(Disp)
					if(sum(DispSuc) > 0){
						Dists <- Dists[order(Dists$DispProb, decreasing = T), ]
						Dists$CumSum <- cumsum(Dists$DispProb)						
						Draw <- runif(1, min = 0, max = sum(Dists$DispProb))
						DrawSuc <- NULL
						DrawSuc <- Dists[Draw < Dists$CumSum & Draw >= c(0, Dists$CumSum[ - nrow(Dists)]), ]
					}
				}
# let one emigrating non-breeder immigrate at drawn GridID if it could immigrate
				if(length(DrawSuc) > 0){
					SJL[SJL$GridID == DrawSuc$GridID, ]$NewNB <- SJL[SJL$GridID == DrawSuc$GridID, ]$NewNB + 1
# remove the non-breeder from the list of emigrants
					SJL[SJL$GridID == CurEmi$GridID, ]$EmiNBToNB <- SJL[SJL$GridID == CurEmi$GridID, ]$EmiNBToNB - 1
				} else {
# let emigrant remain at source group
					SJL[SJL$GridID == CurEmi$GridID, ]$NewNB <- SJL[SJL$GridID == CurEmi$GridID, ]$NewNB + 1
# remove the non-breeder from the list of emigrants
					SJL[SJL$GridID == CurEmi$GridID, ]$EmiNBToNB <- SJL[SJL$GridID == CurEmi$GridID, ]$EmiNBToNB - 1
				}
			}			
		}
	}
	return(SJL)
}

################################################################## Function to identify emigrating non-breeders ####################################################
####################################################################################################################################################################

FunEmiNBToNB <- function(SJL){
	
########### if there are D (non-breeder in September, these used to be C (level in model))
		
	CurNB <- SJL[SJL$D > 0, ]
	if(nrow(CurNB) > 0) {

		GSTAdd <- FunCalcGSTerr(CurNB)

# the covars and standardize

		CurCovarsEmiNBToNB <- data.frame(GSTAdd$GroupSize, stg3 = factor("C", levels = c("C", "D", "X", "Z")))
		colnames(CurCovarsEmiNBToNB)[1] <- "GroupSize"
		CurCovarsEmiNBToNB$GroupSize <- (CurCovarsEmiNBToNB$GroupSize - MEmiNBToNB_GroupSize) / SEmiNBToNB_GroupSize
	
# simulate from model
		EmiData <- model.matrix(~ stg3 + GroupSize + I(GroupSize^2), data = CurCovarsEmiNBToNB)
		LinPredsEmi <- EmiData %*% MEmiNBToNB_Coef
		PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
		PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurNB$D, prob = PredsEmi)
		PredsEmi <- data.frame(CurNB$GridID, PredsEmi)
		colnames(PredsEmi) <- c("GridID", "EmiNBToNB")

# add to new SJLog
		SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] <- PredsEmi$EmiNBToNB
	
# remove emigrants from D log
		SJL$D[match(PredsEmi$GridID, SJL$GridID)] <- SJL$D[match(PredsEmi$GridID, SJL$GridID)] - PredsEmi$EmiNBToNB
	}

########### if there are C (these used to be D)
		
	CurNB <- SJL[SJL$C > 0, ]
	if(nrow(CurNB) > 0) {

		GSTAdd <- FunCalcGSTerr(CurNB)

# the covars and standardize
		CurCovarsEmiNBToNB <- data.frame(GSTAdd$GroupSize, stg3 = factor("D", levels = c("C", "D", "X", "Z")))
		colnames(CurCovarsEmiNBToNB)[1] <- "GroupSize"
		CurCovarsEmiNBToNB$GroupSize <- (CurCovarsEmiNBToNB$GroupSize - MEmiNBToNB_GroupSize) / SEmiNBToNB_GroupSize
	
# simulate from model
		EmiData <- model.matrix(~ stg3 + GroupSize + I(GroupSize^2), data = CurCovarsEmiNBToNB)
		LinPredsEmi <- EmiData %*% MEmiNBToNB_Coef
		PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
		PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurNB$C, prob = PredsEmi)
		PredsEmi <- data.frame(CurNB$GridID, PredsEmi)
		colnames(PredsEmi) <- c("GridID", "EmiNBToNB")

# add to new SJLog
		SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] <- PredsEmi$EmiNBToNB
# remove emigrants from C log
		SJL$C[match(PredsEmi$GridID, SJL$GridID)] <- SJL$C[match(PredsEmi$GridID, SJL$GridID)] - PredsEmi$EmiNBToNB
	}
	

########### if there are XToNB (retained juvenile)
		
	CurNB <- SJL[SJL$XToNB > 0, ]
	if(nrow(CurNB) > 0) {

		GSTAdd <- FunCalcGSTerr(CurNB)

# the covars and standardize
		CurCovarsEmiNBToNB <- data.frame(GSTAdd$GroupSize, stg3 = factor("X", levels = c("C", "D", "X", "Z")))
		colnames(CurCovarsEmiNBToNB)[1] <- "GroupSize"
		CurCovarsEmiNBToNB$GroupSize <- (CurCovarsEmiNBToNB$GroupSize - MEmiNBToNB_GroupSize) / SEmiNBToNB_GroupSize
	
# simulate from model
		EmiData <- model.matrix(~ stg3 + GroupSize + I(GroupSize^2), data = CurCovarsEmiNBToNB)
		LinPredsEmi <- EmiData %*% MEmiNBToNB_Coef
		PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
		PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurNB$XToNB, prob = PredsEmi)
		PredsEmi <- data.frame(CurNB$GridID, PredsEmi)
		colnames(PredsEmi) <- c("GridID", "EmiNBToNB")

# add to new SJLog
		SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] <- SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] + PredsEmi$EmiNBToNB
# remove emigrants from XToNB log
		SJL$XToNB[match(PredsEmi$GridID, SJL$GridID)] <- SJL$XToNB[match(PredsEmi$GridID, SJL$GridID)] - PredsEmi$EmiNBToNB
	}
	

########### if there are ZToNB (dispersed juvenile)
		
	CurNB <- SJL[SJL$ZToNB > 0, ]
	if(nrow(CurNB) > 0) {

		GSTAdd <- FunCalcGSTerr(CurNB)

# the covars and standardize
		CurCovarsEmiNBToNB <- data.frame(GSTAdd$GroupSize, stg3 = factor("Z", levels = c("C", "D", "X", "Z")))
		colnames(CurCovarsEmiNBToNB)[1] <- "GroupSize"
		CurCovarsEmiNBToNB$GroupSize <- (CurCovarsEmiNBToNB$GroupSize - MEmiNBToNB_GroupSize) / SEmiNBToNB_GroupSize
	
# simulate from model
		EmiData <- model.matrix(~ stg3 + GroupSize + I(GroupSize^2), data = CurCovarsEmiNBToNB)
		LinPredsEmi <- EmiData %*% MEmiNBToNB_Coef
		PredsEmi <- exp(LinPredsEmi) / (1 + exp(LinPredsEmi) )
		PredsEmi <- rbinom(n = nrow(PredsEmi), size = CurNB$ZToNB, prob = PredsEmi)
		PredsEmi <- data.frame(CurNB$GridID, PredsEmi)
		colnames(PredsEmi) <- c("GridID", "EmiNBToNB")

# add to new SJLog
		SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] <- SJL$EmiNBToNB[match(PredsEmi$GridID, SJL$GridID)] + PredsEmi$EmiNBToNB
# remove emigrants from ZToNB log
		SJL$ZToNB[match(PredsEmi$GridID, SJL$GridID)] <- SJL$ZToNB[match(PredsEmi$GridID, SJL$GridID)] - PredsEmi$EmiNBToNB
	}
	

	return(SJL)
}
