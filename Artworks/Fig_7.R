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
pred <- function (pars, Species, Dose, idata, route = 'oral', tinterval = 24, Dtimes = 1, endtime = 24*365*5, tt = 24, atol = 1E-10, rtol = atol/2, N = 1000) {
    
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
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, tt) 
    
    ## Simulation
    out <-  mod %>% 
            param (pars) %>%                 
            update(atol = atol, rtol = rtol, maxsteps = 50000) %>%  
            mrgsim (idata = idata, ev = ex, tgrid = tsamp)
    
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    
    if (Species == "BC") {
        outdf = cbind.data.frame (ID    = out$ID,
                                  Time   = out$time, 
                                  CP     = out$Plasma*1000, 
                                  CL     = out$Liver*1000,
                                  CK     = out$Kidney*1000, 
                                  CM     = out$Muscle*1000
                                  )
        
    } else {
        outdf = cbind.data.frame (ID     = out$ID,
                                  Time   = out$time, 
                                  CP     = out$Plasma*1000, 
                                  CL     = out$Liver*1000,
                                  CK     = out$Kidney*1000, 
                                  CM     = out$Muscle*1000,
                                  CMilk  = out$Milk*1000,
                                  Urine  = out$Urine,
                                  AUCCP  = out$AUC_CP
                                  )
  
    }
    

    return (outdf)
    
}



## Calculation of idata
N=100
idata_BC_PFOA   <- MCsim (N=N, pars_BC_PFOA, pars_SD_BC)
idata_BC_PFOS   <- MCsim (N=N, pars_BC_PFOS, pars_SD_BC)
idata_BC_PFHXS  <- MCsim (N=N, pars_BC_PFHXS, pars_SD_BC)


idata_DC_PFOA   <- MCsim (N=N, pars_DC_PFOA, pars_SD_DC)
idata_DC_PFOS   <- MCsim (N=N, pars_DC_PFOS, pars_SD_DC)
idata_DC_PFHXS  <- MCsim (N=N, pars_DC_PFHXS, pars_SD_DC)


## Exposure parameters --------------------------------------------------------
## Please refer to Excel 1 for details about exposure scenarios
## 2-years and 3-years exposure time for cattle and cow, respectively

WDI_fun <- function(Chemical, Species, Dose) {
    
if (Chemical == "PFOA") {
    if (Species == "BC") {
        
        idata = idata_BC_PFOA
        pars = pars_BC_PFOA
        Action_level_M = 10
        Action_level_L = 78
        Action_level_K = 1719
        
    } else {
        idata = idata_DC_PFOA
        pars = pars_DC_PFOA
        Action_level_M = 10
        Action_level_L = 78
        Action_level_K = 1719
        Action_level_Milk = 5.9
        
    }
 }  else if (Chemical == "PFOS") {
        if (Species == "BC") {
            
            idata = idata_BC_PFOS
            pars = pars_BC_PFOS
            Action_level_M = 10
            Action_level_L = 78
            Action_level_K = 1719

        } else {
            idata = idata_DC_PFOS
            pars = pars_DC_PFOS
            Action_level_M = 10
            Action_level_L = 78
            Action_level_K = 1719
            Action_level_Milk = 5.92
        }
        
  } else {
       if (Species == "BC") {
            
            idata = idata_BC_PFHXS
            pars = pars_BC_PFHXS
            Action_level_M = 2
            Action_level_L = 16
            Action_level_K = 343

        } else {
            idata = idata_DC_PFHXS
            pars = pars_DC_PFHXS
            Action_level_M = 2
            Action_level_L = 16
            Action_level_K = 343
            Action_level_Milk = 1.2
        }
    }
    
 Sim <- pred (N=N, pars = pars,  Species = Species, 
              Dose  = Dose, Dtimes  = 365*2, 
              endtime= 365*24*12, idata = idata, tt = 24)
 
 ## Estimated the Withdraw intervals
 WDI <-Sim %>% group_by (Time) %>% summarise (
        M99 = quantile(CM, probs  = 0.99),
        L99 = quantile(CL, probs  = 0.99),
        K99 = quantile(CK, probs  = 0.99),
        Milk99 = ifelse(Species =="DC", quantile(CMilk, probs  = 0.99), quantile(CM, probs  = 0.99))
 )
 
 WDIs_M  <- WDI  %>% filter(Time > 24*12*30*2 & round(M99,3) <= Action_level_M)%>%
            mutate(WDIs = ((Time/24)-12*30*2)+1)%>%select(WDIs)%>%min()
 WDIs_L  <- WDI  %>% filter(Time > 24*12*30*2 & round(L99,3) <= Action_level_L)%>%
            mutate(WDIs = ((Time/24)-12*30*2)+1)%>%select(WDIs)%>%min()
 WDIs_K  <- WDI  %>% filter(Time > 24*12*30*2 & round(K99,3) <= Action_level_K)%>%
            mutate(WDIs = ((Time/24)-12*30*2)+1)%>%select(WDIs)%>%min()
 
 WDIs_Milk  <- ifelse(Species =="DC", WDI  %>% filter(Time > 24*12*30*2 & round(Milk99,3) <= Action_level_Milk)%>%
               mutate(WDIs = ((Time/24)-12*30*2)+1)%>%select(WDIs)%>%min(), 0)
 
 
 return(list(WDIs_M=WDIs_M,
             WDIs_L=WDIs_L,
             WDIs_K=WDIs_K,
             WDIs_Milk = WDIs_Milk, 
             Sim=Sim))
    
}    

## Estimation of WDIS for beef cattle
## High exposure
Dose_BC_PFOA = 4.93E-5
Dose_DC_PFOA = 7.34E-5

Dose_BC_PFOS = 1.93E-4
Dose_DC_PFOS = 2.88E-4
    
Dose_BC_PFHXS = 1.61E-4
Dose_DC_PFHXS = 2.40E-4   
    
    
# PFOA
WDIs_PFOA_BC_H    = WDI_fun(Chemical="PFOA", Species="BC", Dose =Dose_BC_PFOA)
WDIs_BC_PFOA_M_H  = WDIs_PFOA_BC_H[[1]]
WDIs_BC_PFOA_L_H  = WDIs_PFOA_BC_H[[2]]
WDIs_BC_PFOA_K_H  = WDIs_PFOA_BC_H[[3]]
Pdata_BC_PFOA_H   = WDIs_PFOA_BC_H[[5]]

# PFOS
WDIs_PFOS_BC_H    = WDI_fun(Chemical="PFOS", Species="BC", Dose = Dose_BC_PFOS)
WDIs_BC_PFOS_M_H  = WDIs_PFOS_BC_H[[1]]
WDIs_BC_PFOS_L_H  = WDIs_PFOS_BC_H[[2]]
WDIs_BC_PFOS_K_H  = WDIs_PFOS_BC_H[[3]]
Pdata_BC_PFOS_H   = WDIs_PFOS_BC_H[[5]]

# PFHXS
WDIs_PFHXS_BC_H    = WDI_fun(Chemical="PFHXS", Species="BC", Dose = Dose_BC_PFHXS)
WDIs_BC_PFHXS_M_H  = WDIs_PFHXS_BC_H[[1]]
WDIs_BC_PFHXS_L_H  = WDIs_PFHXS_BC_H[[2]]
WDIs_BC_PFHXS_K_H  = WDIs_PFHXS_BC_H[[3]]
Pdata_BC_PFHXS_H   = WDIs_PFHXS_BC_H[[5]]

## Estimation of WDIS for dairy cows
# PFOA
WDIs_PFOA_DC_H    = WDI_fun(Chemical="PFOA", Species="DC", Dose = Dose_DC_PFOA)
WDIs_DC_PFOA_M_H  = WDIs_PFOA_DC_H[[1]]
WDIs_DC_PFOA_L_H  = WDIs_PFOA_DC_H[[2]]
WDIs_DC_PFOA_K_H  = WDIs_PFOA_DC_H[[3]]
WDIs_DC_PFOA_Milk_H  = WDIs_PFOA_DC_H[[4]]
Pdata_DC_PFOA_H    = WDIs_PFOA_DC_H[[5]]

# PFOS
WDIs_PFOS_DC_H    = WDI_fun(Chemical="PFOS", Species="DC", Dose = Dose_DC_PFOS)
WDIs_DC_PFOS_M_H  = WDIs_PFOS_DC_H[[1]]
WDIs_DC_PFOS_L_H  = WDIs_PFOS_DC_H[[2]]
WDIs_DC_PFOS_K_H  = WDIs_PFOS_DC_H[[3]]
WDIs_DC_PFOS_Milk_H  = WDIs_PFOS_DC_H[[4]]
Pdata_DC_PFOS_H    = WDIs_PFOS_DC_H[[5]]

# PFHXS
WDIs_PFHXS_DC_H     = WDI_fun(Chemical="PFHXS", Species="DC", Dose = Dose_DC_PFHXS)
WDIs_DC_PFHXS_M_H   = WDIs_PFHXS_DC_H[[1]]
WDIs_DC_PFHXS_L_H   = WDIs_PFHXS_DC_H[[2]]
WDIs_DC_PFHXS_K_H   = WDIs_PFHXS_DC_H[[3]]
WDIs_DC_PFHXS_Milk_H   = WDIs_PFHXS_DC_H[[4]]
Pdata_DC_PFHXS_H    = WDIs_PFHXS_DC_H[[5]]

#-------------------------------------------------------------------------------
## Middle exposure
Dose_BC_PFOA_M  = 1.35E-6
Dose_DC_PFOA_M  = 2.01E-6

Dose_BC_PFOS_M  = 5.36E-6
Dose_DC_PFOS_M  = 7.98E-6

Dose_BC_PFHXS_M = 1.97E-6
Dose_DC_PFHXS_M = 2.94E-6   

# PFOA
WDIs_PFOA_BC_M    = WDI_fun(Chemical="PFOA", Species="BC", Dose =Dose_BC_PFOA_M)
WDIs_BC_PFOA_M_M  = WDIs_PFOA_BC_M[[1]]
WDIs_BC_PFOA_L_M  = WDIs_PFOA_BC_M[[2]]
WDIs_BC_PFOA_K_M  = WDIs_PFOA_BC_M[[3]]
Pdata_BC_PFOA_M   = WDIs_PFOA_BC_M[[5]]

# PFOS
WDIs_PFOS_BC_M  = WDI_fun(Chemical="PFOS", Species="BC", Dose = Dose_BC_PFOS_M)
WDIs_BC_PFOS_M_M  = WDIs_PFOS_BC_M[[1]]
WDIs_BC_PFOS_L_M  = WDIs_PFOS_BC_M[[2]]
WDIs_BC_PFOS_K_M  = WDIs_PFOS_BC_M[[3]]
Pdata_BC_PFOS_M    = WDIs_PFOS_BC_M[[5]]

# PFHXS
WDIs_PFHXS_BC_M  = WDI_fun(Chemical="PFHXS", Species="BC", Dose = Dose_BC_PFHXS_M)
WDIs_BC_PFHXS_M_M  = WDIs_PFHXS_BC_M[[1]]
WDIs_BC_PFHXS_L_M  = WDIs_PFHXS_BC_M[[2]]
WDIs_BC_PFHXS_K_M  = WDIs_PFHXS_BC_M[[3]]
Pdata_BC_PFHXS_M    = WDIs_PFHXS_BC_M[[5]]

## Estimation of WDIS for dairy cows
# PFOA
WDIs_PFOA_DC_M    = WDI_fun(Chemical="PFOA", Species="DC", Dose = Dose_DC_PFOA_M)
WDIs_DC_PFOA_M_M  = WDIs_PFOA_DC_M[[1]]
WDIs_DC_PFOA_L_M  = WDIs_PFOA_DC_M[[2]]
WDIs_DC_PFOA_K_M  = WDIs_PFOA_DC_M[[3]]
WDIs_DC_PFOA_Milk_M  = WDIs_PFOA_DC_M[[4]]
Pdata_DC_PFOA_M    = WDIs_PFOA_DC_M[[5]]

# PFOS
WDIs_PFOS_DC_M    = WDI_fun(Chemical="PFOS", Species="DC", Dose = Dose_DC_PFOS_M)
WDIs_DC_PFOS_M_M  = WDIs_PFOS_DC_M[[1]]
WDIs_DC_PFOS_L_M  = WDIs_PFOS_DC_M[[2]]
WDIs_DC_PFOS_K_M = WDIs_PFOS_DC_M[[3]]
WDIs_DC_PFOS_Milk_M  = WDIs_PFOS_DC_M[[4]]
Pdata_DC_PFOS_M    = WDIs_PFOS_DC_M[[5]]

# PFHXS
WDIs_PFHXS_DC_M  = WDI_fun(Chemical="PFHXS", Species="DC", Dose = Dose_DC_PFHXS_M)
WDIs_DC_PFHXS_M_M  = WDIs_PFHXS_DC_M[[1]]
WDIs_DC_PFHXS_L_M  = WDIs_PFHXS_DC_M[[2]]
WDIs_DC_PFHXS_K_M  = WDIs_PFHXS_DC_M[[3]]
WDIs_DC_PFHXS_Milk_M  = WDIs_PFHXS_DC_M[[4]]
Pdata_DC_PFHXS_M    = WDIs_PFHXS_DC_M[[5]]


#######
## Low exposure
Dose_BC_PFOA_L  = 1.87E-7
Dose_DC_PFOA_L  = 2.78E-7

Dose_BC_PFOS_L  = 1.65E-7
Dose_DC_PFOS_L  = 2.45E-7

Dose_BC_PFHXS_L = 3.56E-7
Dose_DC_PFHXS_L = 5.30E-7   

# PFOA
WDIs_PFOA_BC_L    = WDI_fun(Chemical="PFOA", Species="BC", Dose =Dose_BC_PFOA_L)
WDIs_BC_PFOA_M_L  = WDIs_PFOA_BC_L[[1]]
WDIs_BC_PFOA_L_L  = WDIs_PFOA_BC_L[[2]]
WDIs_BC_PFOA_K_L  = WDIs_PFOA_BC_L[[3]]
Pdata_BC_PFOA_L   = WDIs_PFOA_BC_L[[5]]

# PFOS
WDIs_PFOS_BC_L  = WDI_fun(Chemical="PFOS", Species="BC", Dose = Dose_BC_PFOS_L)
WDIs_BC_PFOS_M_L  = WDIs_PFOS_BC_L[[1]]
WDIs_BC_PFOS_L_L  = WDIs_PFOS_BC_L[[2]]
WDIs_BC_PFOS_K_L  = WDIs_PFOS_BC_L[[3]]
Pdata_BC_PFOS_L    = WDIs_PFOS_BC_L[[5]]

# PFHXS
WDIs_PFHXS_BC_L  = WDI_fun(Chemical="PFHXS", Species="BC", Dose = Dose_BC_PFHXS_L)
WDIs_BC_PFHXS_M_L  = WDIs_PFHXS_BC_L[[1]]
WDIs_BC_PFHXS_L_L  = WDIs_PFHXS_BC_L[[2]]
WDIs_BC_PFHXS_K_L  = WDIs_PFHXS_BC_L[[3]]
Pdata_BC_PFHXS_L    = WDIs_PFHXS_BC_L[[5]]

## Estimation of WDIS for dairy cows
# PFOA
WDIs_PFOA_DC_L    = WDI_fun(Chemical="PFOA", Species="DC", Dose = Dose_DC_PFOA_L)
WDIs_DC_PFOA_M_L  = WDIs_PFOA_DC_L[[1]]
WDIs_DC_PFOA_L_L  = WDIs_PFOA_DC_L[[2]]
WDIs_DC_PFOA_K_L  = WDIs_PFOA_DC_L[[3]]
WDIs_DC_PFOA_Milk_L  = WDIs_PFOA_DC_L[[4]]
Pdata_DC_PFOA_L    = WDIs_PFOA_DC_L[[5]]

# PFOS
WDIs_PFOS_DC_L    = WDI_fun(Chemical="PFOS", Species="DC", Dose = Dose_DC_PFOS_L)
WDIs_DC_PFOS_M_L  = WDIs_PFOS_DC_L[[1]]
WDIs_DC_PFOS_L_L  = WDIs_PFOS_DC_L[[2]]
WDIs_DC_PFOS_K_L = WDIs_PFOS_DC_L[[3]]
WDIs_DC_PFOS_Milk_L  = WDIs_PFOS_DC_L[[4]]
Pdata_DC_PFOS_L    = WDIs_PFOS_DC_L[[5]]

# PFHXS
WDIs_PFHXS_DC_L  = WDI_fun(Chemical="PFHXS", Species="DC", Dose = Dose_DC_PFHXS_L)
WDIs_DC_PFHXS_M_L  = WDIs_PFHXS_DC_L[[1]]
WDIs_DC_PFHXS_L_L  = WDIs_PFHXS_DC_L[[2]]
WDIs_DC_PFHXS_K_L  = WDIs_PFHXS_DC_L[[3]]
WDIs_DC_PFHXS_Milk_L  = WDIs_PFHXS_DC_L[[4]]
Pdata_DC_PFHXS_L    = WDIs_PFHXS_DC_L[[5]]


## MAKE PLOT FOR Beef cattle and cows--------------------------------------------
## Beef cattle 
## Plot for beef muscle
p_BC_PFOS_L  <- MCplot(data=Pdata_BC_PFOS_L, target = "Muscle", TOL = 10, tdoses =  12*30*2, Y.limit = 1e-8, X.limit = 1881)
p_BC_PFOS_M  <- MCplot(data=Pdata_BC_PFOS_M, target = "Muscle", TOL = 10, tdoses =  12*30*2, Y.limit = 1e-8, X.limit = 1881)
p_BC_PFOS_H  <- MCplot(data=Pdata_BC_PFOS_H, target = "Muscle", TOL = 10, tdoses =  12*30*2, Y.limit = 1e-8, X.limit = 2881)

p_BC_PFOS_L +  p_BC_PFOS_M + p_BC_PFOS_H

## Dairy cows
## Plot for milk
p_DC_PFOS_L <- MCplot(data=Pdata_DC_PFOS_L, Species = 'DC',  target = "Milk",TOL = 5.9, tdoses = 30*12*2, Y.limit = 1e-4, X.limit = 365*3,color = 'bisque3') 
p_DC_PFOS_M <- MCplot(data=Pdata_DC_PFOS_M, Species = 'DC',  target = "Milk",TOL = 5.9, tdoses = 30*12*2, Y.limit = 1e-4, color = 'bisque3') 
p_DC_PFOS_H <- MCplot(data=Pdata_DC_PFOS_H, Species = 'DC',  target = "Milk",TOL = 5.9, tdoses = 30*12*2, Y.limit = 1e-4, color = 'bisque3') 

p_DC_PFOS_L + p_DC_PFOS_M + p_DC_PFOS_H

## Merge plot
(p_BC_PFOS_L +  p_BC_PFOS_M + p_BC_PFOS_H) / (p_DC_PFOS_L + p_DC_PFOS_M+p_DC_PFOS_H)



ggsave("Fig 6.tiff",scale = 1,
       plot = (p_BC_PFOS_L +  p_BC_PFOS_M + p_BC_PFOS_H) / (p_DC_PFOS_L + p_DC_PFOS_M + p_DC_PFOS_H),
       width = 35, height = 15, units = "cm", dpi=320)




