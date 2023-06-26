#------------------------------------------------------------------------------
# - Author        : Wei-Chun Chou
# - Supervisor    : Zhoumeng Lin
# - Date          : June, 2023
# - Abbreviation  : BC, Beef cattle; DC, Dairy cow  
#------------------------------------------------------------------------------



## Translated into mrgsolve based on the code published by Chou and Lin (Environ Int, 2019)
library(mrgsolve)   ## for PBPK modeling
library(ggplot2)    ## for graph demonstration
library(dplyr)      ## for dataframe manipulation
library(FME)        ## for sensitive parameter optimization
library(minpack.lm) ## selecting parameter for optimization, combining with FME
library(patchwork)  ## R-package for combination of seperate ggplot

## Loading the R code
source (file = 'Beef_PBPK.R')
source (file = 'Cow_PBPK.R')
source (file = 'Function.R')      # Loading the parameters

## Input the fit rds file
Fit_BC_PFOS  = readRDS(file = 'Fit_BC_PFOS.rds')
Fit_BC_PFOA  = readRDS(file = 'Fit_BC_PFOA.rds')
Fit_BC_PFHXS = readRDS(file = 'Fit_BC_PFHXS.rds')
Fit_DC_PFOS  = readRDS(file = 'Fit_DC_PFOS.rds')
Fit_DC_PFOA  = readRDS(file = 'Fit_DC_PFOA.rds')
Fit_DC_PFHXS = readRDS(file = 'Fit_DC_PFHXS.rds')

# Extract best-fitted parameters
pars_BC_PFOS <- Fit_BC_PFOS$par
pars_BC_PFOA <- Fit_BC_PFOA$par
pars_BC_PFHXS <- Fit_BC_PFHXS$par

pars_DC_PFOS <- Fit_DC_PFOS$par
pars_DC_PFOA <- Fit_DC_PFOA$par
pars_DC_PFHXS <- Fit_DC_PFHXS$par



## Loading the mrgsolve model
mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)

## Read data
data    <- read.csv(file = "Data.csv")


## Define the prediction function
pred <- function(fixpars = (mod_bc %>% param), pars, BW, tdose, dose, tt=1, Species = "BC", atol = 1E-8, rtol = atol/2, endtime = 24*365*1) {
    
    ## Get out of log domain
    pars <- lapply(pars, exp)           ## Return a list of exp (Parameters for gestational model) from log scale
    mod <- switch(Species,  
                  'BC' = mod_bc,
                  'DC' = mod_dc)
    
    
    ## Exposure scenario for gestational exposure
    BW          = BW                    ## Body weight based on 3-4 yr Belted Galloway cow (https://beltie.org/PDFs/2017/Belted-galloway-Breeders-Manual-Feb-2017.pdf)
    tinterval   = 24                    ## Time interval; 
    TDOSE       = tdose                 ## Total dosing/Dose times; Repeat drinking water dose for 2 yrs
    DOSE        = dose                  ## Input drinking water dose  
    DOSEoral    = DOSE*BW               ## Amount of drinking water dose
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    ex.oral <- ev (ID   = 1,            ## One individual
                   time = 0,            ## Dossed start time 
                   amt  = DOSEoral,     ## Amount of dose 
                   ii   = tinterval,    ## Time interval
                   addl = TDOSE - 1,    ## Additional dosing 
                   cmt  = "AST",        ## The dosing compartment: AST Stomach  
                   replicate = FALSE)   ## No replicate
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, tt) ## Simulation time post 2-yr exposure
    
    ## Simulation of exposure scenaior (Repeated oral dose to 1/2/3/5/10 mg/kg)
    out <- 
        mod %>%  
        param (fixpars) %>%
        param (pars) %>%                 
        update(atol = atol, rtol = rtol, maxsteps = 50000) %>%  
        mrgsim_d (data = ex.oral, tgrid = tsamp) 
    
    
    outdf <- out %>% as.data.frame
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    if (Species=='BC') {
        outdf = cbind.data.frame (Time    = outdf$time, 
                                  CP      = (outdf$Plasma)*1000, 
                                  CL      = (outdf$Liver)*1000,
                                  CK      = (outdf$Kidney)*1000, 
                                  CM      = (outdf$Muscle)*1000,
                                  CU      = outdf$Urine
        )
    
    }else {
        outdf = cbind.data.frame (Time    = outdf$time, 
                                  CP      = outdf$Plasma*1000,
                                  CMilk   = outdf$Milk*1000,
                                  CL      = outdf$Liver*1000,
                                  CK      = outdf$Kidney*1000,
                                  CM      = outdf$Muscle*1000)
    }
    
    return (outdf) 
}


##'''''''''''''''''''' Beef Cattle '''''''''''''''''''''''''''''''''''''''''''''
# PFOS calibration data set  ---------------------------------------------------
# A2: Lupton et al., 2014; single oral exposure at dose of 8 mg/kg; P, L, Mu, K, Sp, Lu
# A3: Lupton et al., 2015; single oral exposure at dose of 0.098 and 9.01 mg/kg; P, L, Mu, K
# A5: Drew, 2019, 2021; feeding with water contaminating 0.00042 mg/kg; P, L, Mu, K 
#-------------------------------------------------------------------------------


BC_PFOS_data <- rbind.data.frame(
    ## Study 1: Guruge et al., 2008
    data_A1_P     <- data %>% filter(Study ==1 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_P", Chemical = 'PFOS'),
    
    ## Study 2: Lupton et al., 2014
    data_A2_P     <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_P", Chemical = 'PFOS'),
    
    data_A2_Mu    <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu")%>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_Mu", Chemical = 'PFOS'),
    
    data_A2_L     <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_L", Chemical = 'PFOS'),
    
    data_A2_K     <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_K", Chemical = 'PFOS'),
    
    data_A2_U     <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "U") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_U", Chemical = 'PFOS'),
    
    ## Study 3: Lupton et al., 2015
    data_A3_P_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "P" & Dose == 0.098) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_P_1", Chemical = 'PFOS'),
    
    data_A3_P_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "P" & Dose == 9.1) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_P_2", Chemical = 'PFOS'),
    
    data_A3_L_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "L" & Dose == 0.098) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_L_1", Chemical = 'PFOS'),
    
    data_A3_L_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "L" & Dose == 9.1) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_L_2", Chemical = 'PFOS'),
    
    data_A3_K_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "K" & Dose == 0.098) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_K_1", Chemical = 'PFOS'),
    
    data_A3_K_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "K" & Dose == 9.1) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_K_2", Chemical = 'PFOS'),
    
    data_A3_Mu_1  <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu" & Dose == 0.098) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_Mu_1", Chemical = 'PFOS'),
    
    data_A3_Mu_2  <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu" & Dose == 9.1) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A3_Mu_2", Chemical = 'PFOS'),
    
    ## Study 5: Drew, 2019, 2021
    data_A5_P     <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_P", Chemical = 'PFOS'),
    
    data_A5_Mu    <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu")%>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_Mu", Chemical = 'PFOS'),
    
    data_A5_L     <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_L", Chemical = 'PFOS'),
    
    data_A5_K     <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_K", Chemical = 'PFOS')
)


# PFOA calibration data set  ---------------------------------------------------
# A4: Lupton et al., 2012; single oral exposure at dose of 8 mg/kg; P, L, Mu, K, Sp, Lu
#-----------------------------------------------------------------------------------
## Study 4: Drew, 2019, 2021

BC_PFOA_data <- rbind.data.frame(
    data_B1_P   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A4_P", Chemical = 'PFOA'),
    data_B1_U   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "U") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A4_U", Chemical = 'PFOA')
)


# PFHxS calibration data set  ---------------------------------------------------
# A4: Lupton et al., 2012; single oral exposure at dose of 8 mg/kg; P, L, Mu, K, Sp, Lu
#-----------------------------------------------------------------------------------

BC_PFHXS_data <- rbind.data.frame(
    
    data_C1_P   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "P" & Time > 15120) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_2_P", Chemical = 'PFHxS'),
    data_C1_L   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_2_L", Chemical = 'PFHxS'),
    data_C1_K   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_2_K", Chemical = 'PFHxS'),
    data_C1_Mu  <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "Mu") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A5_2_Mu", Chemical = 'PFHxS')
    
)

## Merge all BC data

BC_data <- rbind.data.frame(
    BC_PFOS_data,
    BC_PFOA_data,
    BC_PFHXS_data
)




## Define the fix parameters
BC_fixpars_PFOS <- c(
    PL           = 1.03,
    PK           = 0.48,
    PM           = 0.08,
    PRest        = 0.1526,
    Free         = 0.015,
    KabsC        = 2.12,
    KurineC      = 0.01,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif         = 0.001,
    KeffluxC     = 0.1,
    KfecesC      = 0.0137,
    RAF_apical   = 0.001,
    Vmax_apical_invitro = 51803,
    Km_apical    = 64.400,    
    RAF_baso     = 0.005,
    Vmax_baso_invitro = 60.99,
    Km_baso      = 34.500
)

## Define the fixed parameters
BC_fixpars_PFOA <- c(
    PL           = 1.03,
    PK           = 0.48,
    PM           = 0.08,
    PRest        = 0.1526,
    Free         = 0.055,
    KabsC        = 2.12,
    KurineC      = 0.1,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif         = 0.01,
    KeffluxC     = 0.1,
    KfecesC      = 0.0137,
    RAF_apical   = 1e-3,
    Vmax_apical_invitro = 0.90,
    Km_apical    = 0.400,
    RAF_baso     = 0.005,
    Vmax_baso_invitro = 0.9,
    Km_baso      = 34.500
)

BC_fixpars_PFHXS <- c(
    PL           = 0.21,
    PK           = 0.27,
    PM           = 0.062,
    PRest        = 0.715,
    Free         = 0.015,
    KabsC        = 2.12,
    KurineC      = 0.9,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif         = 0.01,
    KeffluxC     = 0.01,
    KfecesC      = 0.0137,
    RAF_apical   = 1e-4,
    Vmax_apical_invitro = 51803,
    Km_apical    = 64.400,    
    RAF_baso     = 1e-3,
    Vmax_baso_invitro = 60.99,
    Km_baso      = 34.500
)


## Prediction-----------------------------------------------------
BC_outdf_A1   <- pred (fixpars=BC_fixpars_PFOS, Species = 'BC', pars = pars_BC_PFOS, BW = 329, dose = 0.00000798, tdose =30*28, tt = 960)
BC_outdf_A2   <- pred (fixpars=BC_fixpars_PFOS, pars = pars_BC_PFOS, BW = 329, dose = 8,       tdose = 1, tt = 1)
BC_outdf_A3_1 <- pred (fixpars=BC_fixpars_PFOS, pars = pars_BC_PFOS, BW = 343, dose = 0.098,   tdose = 1, tt = 1)
BC_outdf_A3_2 <- pred (fixpars=BC_fixpars_PFOS, pars = pars_BC_PFOS, BW = 343, dose = 9.09,    tdose = 1, tt = 1)
BC_outdf_A5   <- pred (fixpars=BC_fixpars_PFOS, pars = pars_BC_PFOS, BW = 364, dose = 0.00042, tdose = 365*1.5, tt =120)
BC_outdf_A4   <- pred (fixpars=BC_fixpars_PFOA, pars = pars_BC_PFOA, BW = 329, dose = 1, tdose = 1, tt = 1, atol = 1E-6, rtol = 1E-3, endtime = 24*30)
BC_outdf_A5_2 <- pred (fixpars=BC_fixpars_PFHXS, pars = pars_BC_PFHXS, BW = 364, dose = 0.00084, tdose = 365*2, tt =120)

BC_pre <- rbind.data.frame(
    
    ## Study 1: Guruge et al., 2008
    Pred_A1_P     <- BC_outdf_A1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "A1_P",  Chemical = 'PFOS'),
    
    ## Study 2: Lupton et al., 2014
    Pred_A2_P     <- BC_outdf_A2  %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A2_P",  Chemical = 'PFOS'),
    Pred_A2_Mu    <- BC_outdf_A2  %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A2_Mu", Chemical = 'PFOS'),
    Pred_A2_L     <- BC_outdf_A2  %>% select (Time = Time, Conc. = CL)  %>% mutate(Study = "A2_L" , Chemical = 'PFOS'),
    Pred_A2_K     <- BC_outdf_A2  %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A2_K",  Chemical = 'PFOS'),
    Pred_A2_U     <- BC_outdf_A2  %>% select (Time = Time, Conc. = CU)  %>% mutate(Study = "A2_U",  Chemical = 'PFOS'),
    
    ## Study 3: Lupton et al., 2015
    Pred_A3_P_1   <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A3_P_1",  Chemical = 'PFOS'),
    Pred_A3_P_2   <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A3_P_2",  Chemical = 'PFOS'),
    Pred_A3_L_1   <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CL)  %>% mutate(Study = "A3_L_1",  Chemical = 'PFOS'),
    Pred_A3_L_2   <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CL)  %>% mutate(Study = "A3_L_2",  Chemical = 'PFOS'),
    Pred_A3_K_1   <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A3_K_1",  Chemical = 'PFOS'),
    Pred_A3_K_2   <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A3_K_2",  Chemical = 'PFOS'),
    Pred_A3_Mu_1  <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A3_Mu_1", Chemical = 'PFOS'),
    Pred_A3_Mu_2  <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A3_Mu_2", Chemical = 'PFOS'),
    
    ## Study 5: Drew, 2019, 2021
    Pred_A5_P     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A5_P", Chemical = 'PFOS'),
    Pred_A5_Mu    <- BC_outdf_A5  %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A5_Mu", Chemical = 'PFOS'),
    Pred_A5_L     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CL)  %>% mutate(Study = "A5_L", Chemical = 'PFOS'),
    Pred_A5_K     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A5_K", Chemical = 'PFOS'),
    
    ## PFOA
    Pred_A4_P     <- BC_outdf_A4  %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A4_P", Chemical = 'PFOA'),
    Pred_A4_U     <- BC_outdf_A4  %>% select (Time = Time, Conc. = CU)  %>% mutate(Study = "A4_U", Chemical = 'PFOA'),
     
    ## PFHxS
    Pred_A5_2_P   <- BC_outdf_A5_2  %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A5_2_P", Chemical = 'PFHxS'),
    Pred_A5_2_Mu  <- BC_outdf_A5_2  %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A5_2_Mu", Chemical = 'PFHxS'),
    Pred_A5_2_L   <- BC_outdf_A5_2  %>% select (Time = Time, Conc. = CL)  %>% mutate(Study = "A5_2_L", Chemical = 'PFHxS'),
    Pred_A5_2_K   <- BC_outdf_A5_2  %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A5_2_K", Chemical = 'PFHxS')
    
) 

# rename the levels
new_levels <- c("Study_1; Plasma",    "Study_2; Kidney",   "Study_2; Liver",    "Study_2; Muscle",   "Study_2; Plasma",
                "Study_2; Urine",     "Study_3_1; Kidney", "Study_3_2; Kidney", "Study_3_1; Liver",  "Study_3_2; Liver",
                "Study_3_1; Muscle",  "Study_3_2; Muscle", "Study_3_1; Plasma", "Study_3_2; Plasma", "Study_4; Plasma",
                "Study_4; Urine",     "Study_5_2; Kidney", "Study_5_2; Liver",  "Study_5_2; Muscle", "Study_5_2; Plasma",
                "Study_5_1; Kidney",  "Study_5_1; Liver",  "Study_5_1; Muscle", "Study_5_1; Plasma")

BC_data$Study <- as.factor(BC_data$Study)
BC_pre$Study  <- as.factor(BC_pre$Study)

levels(BC_data$Study)<- new_levels
levels(BC_pre$Study)<- new_levels




## Make plot
select_list <- c(#"Study_1; Plasma",    "Study_2; Kidney",   "Study_2; Liver",    
                 "Study_2; Muscle",   "Study_2; Plasma",    "Study_2; Urine",     
                 #"Study_3_1; Kidney", "Study_3_2; Kidney", "Study_3_1; Liver",  
                 "Study_3_2; Liver", "Study_3_1; Muscle",  "Study_3_2; Muscle", 
                 "Study_4; Plasma",  #"Study_3_1; Plasma", "Study_3_2; Plasma", 
                 #"Study_4; Urine",     "Study_5_2; Kidney", "Study_5_2; Liver",  
                 "Study_5_2; Muscle", "Study_5_2; Plasma",  #"Study_5_1; Kidney",  
                 "Study_5_1; Liver",  "Study_5_1; Muscle", "Study_5_1; Plasma")



p1 <-
    ggplot(BC_data %>% filter(Study %in% select_list), aes(Time/24, Conc.)) +
    geom_point( colour = "steelblue4", fill = "white", size = 1.5, stroke = 1.2) +
    geom_line  (data = BC_pre%>% filter(Study %in% select_list), aes(Time/24, Conc.), linetype = 2, size = 0.7, colour = "black")  +
    facet_wrap (.~Study, scales = "free",  ncol = 3) +
    labs (x = "", y = "") +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, linewidth=1),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(),
        axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position='none', 
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) 


p1

#save the figure as tiff-format
ggsave("Fig 3-1.tiff",
       plot = p1,
       width = 30, height = 20, units = "cm", dpi=320)



##--------------------------DAIRY COW

DC_PFOS_data <- rbind.data.frame(
    data_A1_P   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_P", Chemical = 'PFOS'),
    
    data_A1_M   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "M") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_M", Chemical = 'PFOS'),
    
    data_A1_L   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_L", Chemical = 'PFOS'),
    data_A1_K   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_K", Chemical = 'PFOS'),
    data_A1_Mu  <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "Mu") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A1_Mu", Chemical = 'PFOS'),
    data_A2_M   <- data %>% filter(Study ==8 & Chem == "PFOS" & Species == 'DC' & Matrix == "M") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "A2_M", Chemical = 'PFOS')
    
)

DC_PFOA_data <- rbind.data.frame(
    data_B1_P   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_P", Chemical = 'PFOA'),
    
    data_B1_M   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "M") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_M", Chemical = 'PFOA'),
    
    data_B1_L   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_L", Chemical = 'PFOA'),
    data_B1_K   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_K", Chemical = 'PFOA'),
    data_B1_Mu  <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "Mu") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_Mu", Chemical = 'PFOA')
    
)

DC_PFHXS_data <- rbind.data.frame(
    data_C1_P   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_P", Chemical = 'PFHXS'),
    
    data_C1_M   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "M") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_M", Chemical = 'PFHXS'),
    
    data_C1_L   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_L", Chemical = 'PFHXS'),
    data_C1_K   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_K", Chemical = 'PFHXS'),
    data_C1_Mu  <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "Mu") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_Mu", Chemical = 'PFHXS')
    
)

## Merge all DC data
DC_data <- rbind.data.frame(
    DC_PFOS_data,
    DC_PFOA_data,
    DC_PFHXS_data
)

## Define the fixed parameters
DC_fixpars_PFOS <- c(
    PL         = 1.6534,
    PK         = 1.35,
    PM        = 0.0772,
    PU        = 0.2,
    PRest      = 0.05,
    Free       = 0.0152,
    KabsC      = 2.12,
    KurineC      = 0.01,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif       = 0.001,
    KeffluxC   = 0.5,
    KfecesC    = 0.002,
    RAF_apical = 0.01,
    Vmax_apical_invitro = 51803,
    Km_apical  = 64.400,    
    RAF_baso   = 0.5,
    Vmax_baso_invitro = 506.96,
    Km_baso    = 34.500,
    PMilkM     = 0.9587,
    Kmilking   =  0.1801   
)

DC_fixpars_PFOA <- c(
    PL           = 1.2,
    PK           = 1.5,
    PM           = 0.0772,
    PU           = 0.2,
    PRest        = 0.05,
    Free         = 0.012,
    KabsC        = 2.12,
    KurineC      = 0.5,
    KbileC       = 0.6,
    KehcC        = 0.002,
    Kdif         = 0.1,
    KeffluxC     = 0.02,
    KfecesC      = 0.06,
    RAF_apical   = 0.0005,
    Vmax_apical_invitro = 0.9470,
    Km_apical    = 52.3,    
    RAF_baso     = 4.07,
    Vmax_baso_invitro = 0.4,
    Km_baso      = 27,
    PMilkM       = 0.9587,
    Kmilking     = 0.5801               
)

DC_fixpars_PFHXS <- c(
    PL         = 0.2,
    PK         = 0.4,
    PM         = 0.0772,
    PU         = 0.2,
    PRest      = 0.5,
    Free       = 0.052,
    KabsC      = 2.12,
    KurineC      = 0.02,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif       = 0.01,
    KeffluxC   = 0.2,
    KfecesC    = 0.01,
    RAF_apical   = 50,
    Vmax_apical_invitro = 0.947,
    Km_apical    = 52.3,    
    RAF_baso     = 4.07,
    Vmax_baso_invitro = 0.04,
    Km_baso      = 27,
    PMilkM       = 0.09587,
    Kmilking     = 0.5801                               
)

## Merge prediction data set
DC_outdf_A1   <- pred (fixpars = DC_fixpars_PFOS, Species='DC', pars = pars_DC_PFOS,  BW = 583, dose = 0.0076, tdose = 28, endtime = 24*50)
DC_outdf_A2   <- pred (fixpars = DC_fixpars_PFOS, Species='DC', pars = pars_DC_PFOS, BW = 587, dose = 0.0004/587, tdose = 545, tt = 1)
DC_outdf_B1   <- pred (fixpars = DC_fixpars_PFOA, Species='DC', pars = pars_DC_PFOA,  BW = 583, dose = 0.002, tdose = 28, endtime = 40*24)
DC_outdf_C1   <- pred (fixpars = DC_fixpars_PFHXS, Species='DC', pars = pars_DC_PFHXS, BW = 583, dose = 0.0046, tdose = 28)

DC_Pre <- rbind.data.frame(
    
    Pred_A1_P     <- DC_outdf_A1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "A1_P",   Chemical = 'PFOS'),
    Pred_A1_M     <- DC_outdf_A1  %>% select(Time = Time, Conc. = CMilk)    %>% mutate(Study = "A1_M",   Chemical = 'PFOS'),
    Pred_A1_L     <- DC_outdf_A1  %>% select(Time = Time, Conc. = CL)   %>% mutate(Study = "A1_L",   Chemical = 'PFOS'),
    Pred_A1_K     <- DC_outdf_A1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study = "A1_K",   Chemical = 'PFOS'),
    Pred_A1_Mu    <- DC_outdf_A1  %>% select(Time = Time, Conc. = CM)  %>% mutate(Study = "A1_Mu",  Chemical = 'PFOS'),
    Pred_A2_M     <- DC_outdf_A2  %>% select(Time = Time, Conc. = CMilk) %>% mutate(Study  = "A2_M",  Chemical = 'PFOS'),
    
    Pred_B1_P     <- DC_outdf_B1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "B1_P",   Chemical = 'PFOA'),
    Pred_B1_M     <- DC_outdf_B1  %>% select(Time = Time, Conc. = CMilk)    %>% mutate(Study = "B1_M",   Chemical = 'PFOA'),
    Pred_B1_L     <- DC_outdf_B1  %>% select(Time = Time, Conc. = CL)   %>% mutate(Study = "B1_L",   Chemical = 'PFOA'),
    Pred_B1_K     <- DC_outdf_B1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study = "B1_K",   Chemical = 'PFOA'),
    Pred_B1_Mu    <- DC_outdf_B1  %>% select(Time = Time, Conc. = CM)  %>% mutate(Study = "B1_Mu",  Chemical = 'PFOA'),
    
    Pred_C1_P     <- DC_outdf_C1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "C1_P",   Chemical = 'PFHxS'),
    Pred_C1_M     <- DC_outdf_C1  %>% select(Time = Time, Conc. = CMilk)    %>% mutate(Study = "C1_M",   Chemical = 'PFHxS'),
    Pred_C1_L     <- DC_outdf_C1  %>% select(Time = Time, Conc. = CL)   %>% mutate(Study = "C1_L",   Chemical = 'PFHxS'),
    Pred_C1_K     <- DC_outdf_C1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study = "C1_K",   Chemical = 'PFHxS'),
    Pred_C1_Mu    <- DC_outdf_C1  %>% select(Time = Time, Conc. = CM)  %>% mutate(Study = "C1_Mu",  Chemical = 'PFHxS')
    
) 


# rename the levels
new_levels <- c("Study_1; Kidney",    "Study_1; Liver",   "Study_1; Milk",    "Study_1; Muscle",   
                "Study_1; Plasma",     "Study_2; Milk",  "Study_3; Kidney", "Study_3; Liver",  
                "Study_3; Milk",    "Study_3; Muscle",  "Study_3; Plasma",  "Study_4; Kidney",
                "Study_4; Liver",    "Study_4; Milk",  "Study_4; Muscle",  "Study_4; Plasma")
              

DC_data$Study <- as.factor(DC_data$Study)
DC_Pre$Study  <- as.factor(DC_Pre$Study)

levels(DC_data$Study)<- new_levels
levels(DC_Pre$Study)<- new_levels


## Make plot
select_list <- c("Study_1; Kidney",   "Study_1; Liver",   "Study_1; Milk",    
                 #"Study_1; Plasma",  "Study_2; Milk",  #"Study_1; Muscle", 
                 "Study_3; Kidney", "Study_3; Liver",  #"Study_3; Milk",    
                 "Study_3; Muscle",  "Study_3; Plasma",  "Study_4; Kidney",
                 "Study_4; Liver",    "Study_4; Milk",  "Study_4; Muscle",  
                 "Study_4; Plasma")

p2 <-
    ggplot(DC_data%>%filter(Study %in% select_list), aes(Time/24, Conc.)) +
    geom_point  (colour = "red", fill = "white", size = 1.5, stroke = 1.2) +
    geom_line  (data = DC_Pre%>%filter(Study %in% select_list), linetype = 2, size = 0.7, colour = "black")  +
    #annotation_logticks() + 
    # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(scale_y_continuous(trans  = scales::log2_trans(),
    #                                       breaks = scales::trans_breaks("log2", function(x) 2^x),
    #                                       labels = scales::trans_format("log2", scales::math_format(2^.x))))+
    #scale_x_continuous(labels = scales::math_format(10^.x)) +
    facet_wrap (.~Study, scales = "free", ncol = 3) + 
    labs (x = "", y = "") +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, linewidth=1),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(),
        axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position='none', 
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) 

p2

#save the figure as tiff-format
ggsave("Fig 4.tiff",
       plot = p2,
       width = 30, height = 20, units = "cm", dpi=320)




