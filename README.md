# R code for an individual-based model of Siberian jay

Code written by Ute Bradter, ute.bradter@nina.no with packages based on PROJ.4. 

For Bradter et al. Habitat suitability models based on opportunistic citizen science data - evaluating forecasts from alternative models versus an individual-based model


##### MARK models required that life-history stages of Siberian jay were represented by a one letter abbreviation and this abbreviation has been used in the IBM:
E: summer breeder observed in March  
F: winter breeder observed in September  
C: summer non-breeder observed in March  
D: winter non-breeder observed in September  
Z: dispersed juvenile observed in September  
X: retained juvenile observed in September  
 
##### Covariates from various sources can have abbreviations that differ from the abbreviation in the associated manuscript:
MeanAge: MAge  
MeanVol: MVol  
PlusDays: PD / MeanDaysAboveZero  
Neighbourhood: NbhAge / PatchVAMeanAge10km  
WinterTemp in HSM: TempJF  
SpringPrec in HSM: PrecAM  

##### Function arguments
YearMin: Start year  
YearMax: End year  
YearsScenarioPeriod: Interval (number of years) at which forest conditions will be updated  
StandsShp: optional shapefile to restrict IBM to a smaller area  
CovarDir: directory from which covariates will be read in  
StartNoRep: first repetition  
NoRep: final repetition  
SJLogInit: initial conditions  
TargetAreaName: Name prefixed to output files  
NamingComment: Name identifying covariates and suffixed to output files  
TempFolder: temporary files from each iteration will be written to this directory  
OutFolder: IBM output will be written to this directory   
