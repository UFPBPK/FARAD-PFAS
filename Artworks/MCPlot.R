#------------------------------------------------------------------------------
# The code is to reproduce the figure 5-6
# Note: Prior to running this code, the code for figure 3 should be run
#------------------------------------------------------------------------------

## Loading required R packages
library(EnvStats)
library(readxl)     ## R-package for importing excel file to R
library(FME)        ## R-package for the function of 'modFit' and 'modCost'
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)
library(ggprism)    ## R-package for theme of 'GraphPad Prism'
library(patchwork)  ## R-package for combination of seperate ggplot

## Input the PBPK model code and cirplot code
source (file = 'Beef_PBPK.R') # Loading the code for beef pbpk model
source (file = 'Cow_PBPK.R')  # Loading the code for cow pbpk model
source (file = 'Function.R')      # Loading the parameters


## Loading the mrgsolve model
mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)

### Define the prediction function
pred <- function (pars, Species, Dose, idata, route = 'oral', tinterval = 24, Dtimes = 1, tt = 365, N = 1000) {
    
    ## Exposure scenarios
    if (is.null(idata$BW) == TRUE) {BW = rep(pars["BW"], N)} else {BW =idata$BW}
    
    
    tinterval   = tinterval
    TDOSE       = Dtimes
    DOSE        = Dose*BW 
    mod         = switch(Species, "BC" = mod_bc, "DC" = mod_dc)
  
        
    if (route == 'oral') {
        ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii = tinterval,
                    addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
        ex <- ev_1
        
    } else {
        ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii = tinterval,
                    addl = TDOSE - 1, cmt  = "APL", replicate = FALSE)
        ex <- ev_1
    }
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*tt, 24*30) 
    
    ## Simulation
    out <-  mod %>% 
            param (pars) %>%                 
            update(atol = 1E-3, rtol = 1E-3, maxsteps = 50000) %>%  
            mrgsim (idata = idata, ev = ex, tgrid = tsamp)
    
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    
    if (Species == "BC") {
        outdf = cbind.data.frame (ID    = out$ID,
                                  Time   = out$time, 
                                  CP     = out$Plasma, 
                                  CL     = out$Liver,
                                  CK     = out$Kidney, 
                                  CM     = out$Muscle,
                                  Urine  = out$Urine,
                                  AUCCP  = out$AUC_CP,
                                  AUCCM  = out$AUC_CM,
                                  AUCCL  = out$AUC_CL,
                                  AUCCK  = out$AUC_CK
        )
        
    } else {
        outdf = cbind.data.frame (ID     = out$ID,
                                  Time   = out$time, 
                                  CP     = out$Plasma, 
                                  CL     = out$Liver,
                                  CK     = out$Kidney, 
                                  CM     = out$Muscle,
                                  CMilk  = out$Milk,
                                  Urine  = out$Urine,
                                  AUCCP  = out$AUC_CP,
                                  AUCCM  = out$AUC_CM,
                                  AUCCL  = out$AUC_CL,
                                  AUCCK  = out$AUC_CK) 
    }
    

    return (outdf)
    
}


#
## Calculation of idata
idata_BC_PFOS <- MCsim (pars_BC_PFOS, pars_SD_BC)
idata_DC_PFOS <- MCsim (pars_DC_PFOS, pars_SD_DC)


## Exposure parameters --------------------------------------------------------
# Cow's parameters 
# Ref.1: Maine CDC (https://www.maine.gov/dep/spills/topics/pfas/Agronomic-Pathway-Soil-Screening-Levels-Soil-Fodder-Cows-Milk-09.16.20.pdf)
# Ref.2: EPA PRG (https://epa-prgs.ornl.gov/radionuclides/)
# Fodder intake rate: 13.2 kg/day in hay; 4.1 kg/day in corn silage (ref. 1)
# Water intake rate: 92 L/day (EPA PRG default value from ref. 1 and ref. 2

# Cattle's parameters
# Ref. 3: Beef Cattle Water Requirements and Source Management (https://reurl.cc/Zr2DRl)
# Fodder intake: 
# Water intake rate: 21.9554/day (calculated from 5.8 gallons of water per day)

## Cattle
Sim_PFOS_BC <-pred (pars = pars_BC_PFOS, Species = 'BC', Dose = 0.00042, Dtimes  = 6*30, tt= 365*2, idata = idata_BC_PFOS)
Sim_PFOS_DC <-pred (pars = pars_DC_PFOS, Species = 'DC', Dose = 0.00042 , Dtimes = 6*30, tt= 365*2, idata = idata_DC_PFOS)

## Plot
p_PG_C_A1_L<- MCplot(data=Sim_PFOS_BC, target = "Muscle", TOL = 0.0073, tdoses = 6*30, Y.limit = 1e-8)

p_PG_C_A1_L 

p_PG_C_A1_K<- MCplot(data=Sim_PFOS_DC, Species = 'DC', target = "Milk",TOL = 0.002, tdoses = 365, Y.limit = 1e-8, color = 'bisque3') 

p_PG_C_A1_K




