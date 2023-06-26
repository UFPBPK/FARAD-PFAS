
## MODEL EVALUATION (PFOA and PFOS for beef cattle)-----------------------------



## MODEL SIMULATION ------------------------------------------------------------
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


# Define the mrg prediction function
### Define the prediction function
pred <- function (pars, Species, Dose, idata, route = 'oral', tinterval = 24, Dtimes = 1, tt = 24, N = 1000) {
    
    ## Exposure scenarios
    BW = rep(512, N)
    
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
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*tt, 24) 
    
    ## Simulation
    out <-  mod %>% 
        param (pars) %>%                 
        update(atol = 1E-6, rtol = 1E-3, maxsteps = 50000) %>%  
        mrgsim (idata = idata, ev = ex, tgrid = tsamp)
    
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    
    if (Species == "BC") {
        PDat <- out %>% group_by(time)%>%
                summarise( M_Plasma  = median(Plasma),
                           SD_Plasma = sd(Plasma),
                           M_Liver   = median(Liver),
                           SD_Liver  = sd(Liver),
                           M_Muscle  = median(Muscle),
                           SD_Muscle = sd(Muscle),
                           M_Kidney  = median(Kidney),
                           SD_Kidney = sd(Kidney)
                )
        
    } else {
        
        PDat <- out %>% group_by(time)%>%
            summarise( M_Plasma  = median(Plasma),
                       SD_Plasma = sd(Plasma),
                       M_Liver   = median(Liver),
                       SD_Liver  = sd(Liver),
                       M_Muscle  = median(Muscle),
                       SD_Muscle = sd(Muscle),
                       M_Kidney  = median(Kidney),
                       SD_Kidney = sd(Kidney),
                       M_Milk    = median(Milk),
                       SD_Milk   = sd(Milk),
                       
            )
    }
    
    return (PDat)
    
}


## 
obs.muscle <- cbind.data.frame (
    PFOA  = c(0.011, 0.00007),
    PFOS  = c(0.008, 0.00003),
    unit =  "ng/g"
)

obs.lver <- cbind.data.frame (
    PFOA  = c(0.264, 0.0089),
    PFOS  = c(0.617, 0.051),
    unit =  "ng/g"
)

obs.milk <- cbind.data.frame (
    PFOA  = c(0.016, NULL),
    PFOS  = c(0.024, NULL),
    unit =  "ng/g"
)


rownames(obs.muscle) = rownames(obs.lver) = c("Mean", "SE") 
rownames(obs.milk)  = c("Mean")

## Exposure parameters /////////////////////////////////////////////////////////
# Cow's parameters-------------- 
# Fodder intake rate: 13.2 kg/day in hay; 4.1 kg/day in corn silage (ref. 1)
# Water intake rate: 92 L/day (EPA PRG default value from ref. 1 and ref. 2

# Cattle's parameters-----------
# Fodder intake: 
# Water intake rate: 21.9554/day (calculated from 5.8 gallons of water per day)
# Water intake rate: 0.14 L/kg/d (NRC, 1996)


# PFAS levels in literature-------
# Ref.4 detected PFAS levels in 293 sames of muscle and beef in Xinjiang, China (Table 1); doi: 10.3390/ijerph14090970
# Ref.5 detected PFAS levels in milk in Xinjiang, china (Table 2); doi: 10.3390/ijerph13101037

# - Water samples in Xinjiang area (ref.6: Table S5): PFOA (0.75 ng/L), PFOS (0.34 ng/L), PFHXS (0.25) 
# - Soil samples: 

# Calculation----------------
# Cattle:
# PFOA exposure from water consumption = 0.14 (L/kg/d) X 0.75 (ng/L) = 0.105 (ng/kg/d)
# PFOS exposure from water consumption = 0.14 (L/kg/d) X 0.34 (ng/L) = 0.0476 (ng/kg/d)

# Cows:
# PFOA exposure from soil consumption = 0.5 (kg/day) x 2.3e-6 (mg/kg) = 1.15e-06 (mg/kg)
# PFOS exposure from soil consumption = 0.5 (kg/day) x 0.209e-6 (mg/kg) = 1.045e-07 (mg/kg)

# Reference list-------------
# Ref.1: Maine CDC (https://www.maine.gov/dep/spills/topics/pfas/Agronomic-Pathway-Soil-Screening-Levels-Soil-Fodder-Cows-Milk-09.16.20.pdf)
# Ref.2: EPA PRG (https://epa-prgs.ornl.gov/radionuclides/)
# Ref.3: Beef Cattle Water Requirements and Source Management (https://reurl.cc/Zr2DRl)
# Ref.4: Wang et al. (2017); doi: 10.3390/ijerph14090970
# Ref.5: Xing et al. (2016); doi: 10.3390/ijerph13101037
# Ref.6: Li et al. (2019): doi: https://doi.org/10.1016/j.envint.2018.11.036


#//////////////////////////////////////////////////////////////////////////////

## Calculation of idata
idata_BC_PFOS <- MCsim (pars_BC_PFOS, pars_SD_BC)
idata_BC_PFOA <- MCsim (pars_BC_PFOA, pars_SD_BC)

idata_DC_PFOS <- MCsim (pars_DC_PFOS, pars_SD_DC)
idata_DC_PFOA <- MCsim (pars_DC_PFOA, pars_SD_DC)


## Cattle
Sim_PFOS_BC <-pred (pars = pars_BC_PFOS, Species = 'BC', Dose = (0.0476e-6),  Dtimes = 30*28, tt= 1,  idata = idata_BC_PFOS)
Sim_PFOA_BC <-pred (pars = pars_BC_PFOA, Species = 'BC', Dose = (1.305e-5),   Dtimes = 30*28, tt= 1,  idata = idata_BC_PFOA)
Sim_PFOS_DC <-pred (pars = pars_DC_PFOS, Species = 'DC', Dose = (0.0476e-6),  Dtimes = 30*28, tt= 1,  idata = idata_BC_PFOS)
Sim_PFOA_DC <-pred (pars = pars_DC_PFOA, Species = 'DC', Dose = (1.305e-5),   Dtimes = 30*28, tt= 1,  idata = idata_BC_PFOA)

data_pfos_bc <- Sim_PFOS_BC %>% filter(time == max(time))
data_pfoa_bc <- Sim_PFOA_BC %>% filter(time == max(time))
data_pfos_dc <- Sim_PFOS_DC %>% filter(time == max(time))
data_pfoa_dc <- Sim_PFOA_DC %>% filter(time == max(time))






Pdata_bc <- cbind.data.frame(Mean  = c(data_pfos_bc$M_Muscle*1000, data_pfos_bc$M_Liver*1000,
                                       data_pfoa_bc$M_Muscle*1000, data_pfoa_bc$M_Liver*1000,
                                       data_pfoa_dc$M_Milk*1000, data_pfos_dc$M_Milk*1000,
                                          0.008, 0.617, 0.011, 0.264, 16.2e-3, 24.5e-3),
                               SD    = c(data_pfos_bc$SD_Muscle*1000, data_pfos_bc$SD_Liver*1000,
                                         data_pfoa_bc$SD_Muscle*1000, data_pfoa_bc$SD_Liver*1000,
                                         data_pfoa_dc$SD_Milk*1000, data_pfoa_dc$SD_Milk*1000,
                                          0.0042, 0.27, 0.00098, 0.1246, 0, 0),
                               Type     = c(rep("pred",6), rep("obs",6)),
                               Target   = c(rep(c("PFOS-Muscle", "PFOS-Liver", "PFOA-Muscle", "PFOA-Liver", "PFOA-Milk", "PFOS-Milk"),2)),
                               Chemical = rep(c("PFOS", "PFOS", "PFOA", "PFOA", "PFOA", "PFOS"),2),
                               Species  = c(rep('BC',4), rep('DC',2),rep('BC',4), rep('DC',2))
               
)
                 
Pdata_bc             
       
# Change the colors manually
p<-ggplot(Pdata_bc, aes(Target , Mean, fill=Type)) + 
    geom_bar(stat="identity", color="black", position=position_dodge())+
    geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width=.2,
                  position=position_dodge(.9)) +

    theme_minimal()
    

# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))







