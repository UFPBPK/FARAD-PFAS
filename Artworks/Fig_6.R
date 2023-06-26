## ----------------------------------------------------------------------------
# - Author        : Wei-Chun Chou
# - Supervisor    : Zhoumeng Lin
# - Date          : Started in February 2022, finalized in June 2023
# - Abbreviation  : BC, Beef cattle; DC, Dairy cow  
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The code is to run local sensitivity analysis
#------------------------------------------------------------------------------

## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  
library(ggplot2)
library(tidyr)
library(patchwork)

## Input the PBPK model code and cirplot code
source (file = 'Beef_PBPK.R') # Loading the code for beef pbpk model
source (file = 'Cow_PBPK.R')  # Loading the code for cow pbpk model
source (file = 'Function.R')      # Loading the parameters


## Loading the mrgsolve model
mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)

### Define the prediction function
pred <- function(Species, pars, tdose, dose, tt=1, atol = 1E-6, rtol = atol/2) {
    
    ## Get out of log domain
    pars <- lapply(pars, exp)           ## Return a list of exp (Parameters for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    BW          = if_else(Species == 'BC', 512, 500)                    ## Body weight; 
    tinterval   = 24                    ## Time interval; 
    TDOSE       = tdose                 ## Total dosing/Dose times; 
    DOSE        = dose                  ## Input drinking water dose  
    DOSEoral    = DOSE*BW              ## Amount of drinking water dose
    mod         = switch(Species, "BC" = mod_bc, "DC" = mod_dc)
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    ex.oral <- ev (ID   = 1,            ## One individual
                   amt  = DOSEoral,     ## Amount of dose 
                   ii   = tinterval,    ## Time interval
                   addl = TDOSE - 1,    ## Additional dosing 
                   cmt  = "AST",        ## The dosing compartment: AST Stomach  
                   replicate = FALSE)   ## No replicate
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 4*116, tt) ## Simulation time post 2-yr exposure
    
    ## Simulation of exposure scenario 
    out <- 
        mod %>% 
        param (pars) %>%                 
        update(atol = atol, rtol = rtol, maxsteps = 50000) %>%  
        mrgsim_d (data = ex.oral, tgrid = tsamp) 
    
    

    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    
    if (Species == "BC") {
        outdf = cbind.data.frame (Time   = out$time, 
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
        outdf = cbind.data.frame (Time   = out$time, 
                                  CP     = out$Plasma, 
                                  CL     = out$Liver,
                                  CK     = out$Kidney, 
                                  CM     = out$Muscle,
                                  CMilk  = out$Milk,
                                  Urine  = out$Urine,
                                  AUCCP  = out$AUC_CP,
                                  AUCCM  = out$AUC_CM,
                                  AUCCL  = out$AUC_CL,
                                  AUCCK  = out$AUC_CK,
                                  AUCCMilk = out$AUC_CMilk) 
        }
    
    return (outdf) 
 
}

## Define the sensitivity function
NSC_func <- function (Species, pars, dose, tdose, AUC="AUCCM") {
    # 
  
    n <- length(pars)
    NSC_AUC     <- matrix(NA, nrow = length(pars) , ncol = 1)
    
    for (i in 1:n) {
            pars.new     <- pars %>% replace(i, log(exp((pars[i]))*1.01))
            new          <- pred(Species = Species, pars = pars.new, tdose = tdose, dose = dose)
            ori          <- pred(Species = Species, pars = pars, tdose = tdose, dose = dose)
            delta.pars   <- exp(pars[i])/(exp(pars[i])*0.01)
     
        
        ## Estimated the AUC
        AUC.new   =  new[new$Time==max(new$Time),]-new[new$Time==tdose*24,]
        AUC.ori   =  ori[ori$Time==max(ori$Time),]-ori[ori$Time==tdose*24,]
        
        delta.AUC   =  (AUC.new - AUC.ori)%>%select(AUC)
        AUC.ori_sel =  AUC.ori%>%select(AUC)
        
        NSC_AUC [i, 1]   <- as.numeric((delta.AUC/AUC.ori_sel) * delta.pars)
    }
    
    return (list(NSC_AUC = NSC_AUC))
}


## The exposure scenario for sensitivity analysis: two-year exposure at the daily dose of 0.0001 mg/kg 
## Calculation of NSC for PFOS in Beef Cattle

NSC_BC <- function (AUC) {
    
return( 
    rbind.data.frame(
        NSC_BC_PFOS  <- data.frame(NSC_func (Species = "BC", pars = log(pars_BC_PFOS), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC)) %>% 
                        mutate(Pars = names(pars_BC_PFOS), Chem = "PFOS", Spe = "BC"),
        
        NSC_BC_PFOA  <- data.frame(NSC_func (Species = "BC", pars = log(pars_BC_PFOA), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC)) %>% 
                        mutate(Pars = names(pars_BC_PFOA), Chem = "PFOA", Spe = "BC"),
        
        NSC_BC_PFHXS <- data.frame(NSC_func (Species = "BC", pars = log(pars_BC_PFHXS), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC))%>% 
                        mutate(Pars = names(pars_BC_PFHXS), Chem = "PFHXS", Spe = "BC"))
)
    
}

NSC_BC_AUCCM <- unique(NSC_BC("AUCCM") %>% 
                select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))

NSC_BC_AUCCL <- unique(NSC_BC("AUCCL") %>% 
                select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))

NSC_BC_AUCCK <- unique(NSC_BC("AUCCK") %>% 
                select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))


write.csv(NSC_BC_AUCCM, file = 'a.csv')
write.csv(NSC_DC_AUCCM, file = 'b.csv')
write.csv(NSC_DC_AUCCL, file = 'c.csv')
write.csv(NSC_DC_AUCCK, file = 'd.csv')



## Calculation of NSC for PFOS in Dariy cow

NSC_DC  <- function (AUC) { 
    
return (    
    rbind.data.frame(
        NSC_DC_PFOS  <- data.frame(NSC_func (Species = "DC", pars = log(pars_DC_PFOS), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC)) %>% 
                        mutate(Pars = names(pars_DC_PFOS), Chem = "PFOS", Spe = "DC"),
    
        NSC_DC_PFOA  <- data.frame(NSC_func (Species = "DC", pars = log(pars_DC_PFOA), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC)) %>%
                        mutate(Pars = names(pars_DC_PFOA), Chem = "PFOA", Spe = "DC"),
    
        NSC_DC_PFHXS <- data.frame(NSC_func (Species = "DC", pars = log(pars_DC_PFHXS), 
                        dose = 0.0001, tdose = 365*2, AUC=AUC))%>%
                        mutate(Pars = names(pars_DC_PFHXS), Chem = "PFHXS", Spe = "DC"))
)
}

NSC_DC_AUCCMilk <- unique(NSC_DC("AUCCMilk") %>% 
                   select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))

NSC_DC_AUCCM    <- unique(NSC_DC("AUCCM") %>% 
                   select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))

NSC_DC_AUCCL    <- unique(NSC_DC("AUCCL") %>% 
                   select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))

NSC_DC_AUCCK    <- unique(NSC_DC("AUCCK") %>% 
                   select(Pars, NSC_AUC, Chem, Spe)%>% filter (abs(NSC_AUC) >= 0.1))


## Plot the circle plot
Pdata_BC <- NSC_BC_AUCCM  %>% filter (abs(NSC_AUC) >= 0.1) %>% mutate(value = abs(NSC_AUC)*100) %>% 
            select(Pars, Chem, value)

Pdata_DC <- NSC_DC_AUCCMilk  %>% filter (abs(NSC_AUC) >= 0.1) %>% mutate(value = abs(NSC_AUC)*100) %>% 
            select(Pars, Chem, value)


Pdata_BC$Pars[c(2,3,4,12,13,22:25)] <- c("VLC", "VMC", "VRestC", "VLC", "VRestC","VLC","VKC","VMC","VRestC")
Pdata_DC$Pars[c(2:4, 13:15, 24:26)] <- c("VLC", "VMC", "VRestC", "VLC", "VMC", "VRestC", "VLC" ,"VKC","VRestC")

p1 <- Cirplot (Pdata_BC) + scale_fill_brewer(palette = "Dark2") #scale_fill_manual(values=c("red", "blue", "green"))
p2 <- Cirplot (Pdata_DC) + scale_fill_brewer(palette = "Dark2") #scale_fill_manual(values=c("red", "blue", "green"))

p1 + p2

## Save the figure
ggsave("Fig 5.tiff",scale = 1.1,
       plot = p1+p2,
       width = 35, height = 25, units = "cm", dpi=320)






