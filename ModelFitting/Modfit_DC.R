## Translated into mrgsolve based on the code published by Chou and Lin (Environ Int, 2019)
library(mrgsolve)   ## for PBPK modeling
library(ggplot2)    ## for graph demonstration
library(dplyr)      ## for dataframe manipulation
library(FME)        ## for sensitive parameter optimization
library(minpack.lm) ## selecting parameter for optimization, combining with FME

## Loading the R code
#source (file = 'Beef_PBPK.R')
source (file = 'Cow_PBPK.R')

## Loading the mrgsolve model
#mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)

## Read data
data    <- read.csv(file = "Data.csv")


# PFOS calibration data set  ---------------------------------------------------
# A7: Kowalczyk et al., 2013; L, K, M, Mu, P, U


## Study 7: Kowalczyk et al., 2013
data_A1_P   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.)
data_A1_M   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "M") %>% 
               select (Time = Time, CMilk = Conc.)
data_A1_L   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "L") %>% 
               select (Time = Time, CL = Conc.)
data_A1_K   <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "K") %>% 
               select (Time = Time, CK = Conc.)
data_A1_Mu  <- data %>% filter(Study ==7 & Chem == "PFOS" & Species == 'DC' & Matrix == "Mu") %>% 
               select (Time = Time, CM = Conc.)

## Study 8: Vestergren et al. 2013---------------------------------------------

data_A2_M   <- data %>% filter(Study ==8 & Chem == "PFOS" & Species == 'DC' & Matrix == "M") %>% 
               select (Time = Time, CMilk = Conc.)



## Exposure sceinario 
pred <- function(fixpars = (mod_bc %>% param), pars, BW, tdose, dose, tt=1, atol = 1E-6, rtol = atol/2, endtime = 24*188) {
    
  ## Get out of log domain
  pars <- lapply(pars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  BW          = BW                    ## Body weight based on 3-4 yr Belted Galloway cow (https://beltie.org/PDFs/2017/Belted-galloway-Breeders-Manual-Feb-2017.pdf)
  tinterval   = 24                      ## Time interval; 
  TDOSE       = tdose                      ## Total dosing/Dose times; Repeat oral dose from GD0 - GD40 (weeks)
  DOSE        = dose                     ## Input oral dose (ug/kg/d)
  DOSEoral    = DOSE*BW                 ## Amount of oral dose
  
  # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
  ex.oral <- ev (ID   = 1,              ## One individual
                 time = 0,            ## Dose start time 
                 amt  = DOSEoral,     ## Amount of dose 
                 ii   = tinterval,    ## Time interval
                 addl = TDOSE - 1,    ## Addtional dosing 
                 cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                 replicate = FALSE)   ## No replicate
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, tt) ## Simulation time from GD0 to GD 40
  
  ## Simulation of exposure scenaior (Repeated oral dose to 1/2/3/5/10 mg/kg)
  out <- 
    mod_dc %>% ## Gestational PBPK model
    param (fixpars) %>%
    param (pars) %>%                 ## Update the parameter list with Gpars
    update(atol = atol, rtol=rtol, maxsteps = 50000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
    mrgsim_d (data = ex.oral, tgrid = tsamp) 
  
  
  outdf <- out %>% as.data.frame
  
  ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
  outdf = cbind.data.frame (Time   = outdf$time, 
                            CP = outdf$Plasma*1000,
                            CMilk   = outdf$Milk*1000,
                            CL  = outdf$Liver*1000,
                            CK = outdf$Kidney*1000,
                            CM = outdf$Muscle*1000
  )
  
  return (outdf) # Return Goutdf
}

## Define the fixed parameters
fixpars <- c(
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

## Create a cost function and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
MCost<-function (pars){
  
  outdf_A1   <- pred (fixpars=fixpars, pars = pars, BW = 583, dose = 0.0076, tdose = 28)
  outdf_A2   <- pred (fixpars=fixpars, pars = pars, BW = 587, dose = 0.0004/587, tdose = 545, tt = 1)
  
  cost<- modCost  (model = outdf_A1, obs = data_A1_P, x ="Time", weight = "mean")
  cost<- modCost  (model = outdf_A1, obs = data_A1_L, x ="Time", cost = cost, weight = "mean")
  cost<- modCost  (model = outdf_A1, obs = data_A1_K, x ="Time", cost = cost, weight = "mean")
  cost<- modCost  (model = outdf_A1, obs = data_A1_M, x ="Time", cost = cost, weight = "mean")
  cost<- modCost  (model = outdf_A1, obs = data_A1_Mu, x ="Time", cost = cost, weight = "mean")
  cost<- modCost  (model = outdf_A2, obs = data_A2_M,  x ="Time", cost = cost, weight = "mean")

  
  return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## Initial and/or optimized parameters
theta_init <- log(c(
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
  RAF_apical = 0.001,
  Vmax_apical_invitro = 51803,
  Km_apical  = 64.400,    
  RAF_baso   = 0.5,
  Vmax_baso_invitro = 506.96,
  Km_baso    = 34.500,
  PMilkM     = 0.9587,
  Kmilking   =  0.1801   
))

MCost(theta_init)

## Senstivity function (FME) 
## Check the senstive parameters in the model
# SnsPlasma <- sensFun(func = MCost, parms = theta_init, varscale = 1)
# Sen = summary(SnsPlasma)
# plot(summary(SnsPlasma))

## Selected sensitive parameters;
theta <- theta_init[abs(Sen$Mean) > 1.5*mean(abs(Sen$Mean))]
theta


##
theta <- log(c(
    #PL         = 1.6534,
    PK         = 1.35,
    PM        = 0.0772,
    #PU        = 0.2,
    PRest      = 0.05,
    #Free       = 0.0152,
    #KabsC      = 2.12,
    #KurineC      = 0.01,
    KbileC       = 0.01,
    KehcC        = 0.01,
    Kdif       = 0.001,
    #KeffluxC   = 0.5,
    #KfecesC    = 0.002,
    #RAF_apical = 0.001,
    #Vmax_apical_invitro = 51803,
    Km_apical  = 64.400,   
    #RAF_baso   = 0.5,
    Vmax_baso_invitro = 506.96
    #Km_baso    = 34.500
    #PMilkM     = 0.9587,
    #Kmilking   =  0.1801     
))


## PBPK model fitting 
## Least squares fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
Fit_PFOS <- modFit(f = MCost, p = theta, method ="Marq",
              control = nls.lm.control(nprint = 1))

summary(Fit_PFOS)                                  ## Summary of fit 
exp(Fit_PFOS$par)                                  ## Get the arithmetic value out of the log domain

MCost(Fit_PFOS$par)


##
DC_PFOS<-MCost(Fit_PFOS$par)


PDat <- cbind.data.frame (OBS = DC_PFOS$residuals$obs,
                          PRE = DC_PFOS$residuals$mod,
                          RES = DC_PFOS$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "BC")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat)
summary(fit)

PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
  ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
  geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
  geom_abline (intercept = 0,
               slope     = 1,
               color     ="steelblue4", size = 1, alpha = 0.8) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,4), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,4),labels = scales::math_format(10^.x))


p1

## Time vs. concentrations plot ------------------------------------------------
new_pars <- Fit_PFOS$par
DC_outdf_A1   <- pred (fixpars=fixpars, pars = new_pars, BW = 583, dose = 0.0076, tdose = 28, endtime = 24*50,tt = 1)
DC_outdf_A2   <- pred (fixpars=fixpars, pars = new_pars, BW = 587, dose = 0.0004/587, tdose = 545, tt = 1)

DC_PFOS_pre <- rbind.data.frame(
    
    Pred_A1_P  <- DC_outdf_A1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study    = "A1_P",  Chemical = 'PFOS'),
    Pred_A1_L  <- DC_outdf_A1  %>% select(Time = Time, Conc. = CL)  %>% mutate(Study    = "A1_L",  Chemical = 'PFOS'),
    Pred_A1_K  <- DC_outdf_A1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study    = "A1_K",  Chemical = 'PFOS'),
    Pred_A1_M  <- DC_outdf_A1  %>% select(Time = Time, Conc. = CMilk) %>% mutate(Study  = "A1_M",  Chemical = 'PFOS'),
    Pred_A1_MU <- DC_outdf_A1  %>% select(Time = Time, Conc. = CM) %>% mutate(Study    = "A1_Mu",  Chemical = 'PFOS'),
    Pred_A2_M  <- DC_outdf_A2  %>% select(Time = Time, Conc. = CMilk) %>% mutate(Study  = "A2_M",  Chemical = 'PFOS')
    
)    

DC_OBS_data <- rbind.data.frame(
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



p2 <-
    ggplot(DC_OBS_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = DC_PFOS_pre, aes(Time/24, Conc.))  +
    #annotation_logticks() + 
    # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(scale_y_continuous(trans  = scales::log2_trans(),
    #                                       breaks = scales::trans_breaks("log2", function(x) 2^x),
    #                                       labels = scales::trans_format("log2", scales::math_format(2^.x))))+
    #scale_x_continuous(labels = scales::math_format(10^.x)) +
    facet_wrap (.~Study, scales = "free")


p2





# PFOA calibration data set  ---------------------------------------------------
# A7: Kowalczyk et al., 2013; Repeated exposure at dose of 0.002 mg/kg for 28 days; P, M, Mu, l, k
#----------------------------------------------------------------------------------


## Study 7: Kowalczyk et al., 2013
data_B1_P   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.)
data_B1_M   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "M") %>% 
               select (Time = Time, CMilk = Conc.)
data_B1_L   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "L") %>% 
               select (Time = Time, CL = Conc.)
data_B1_K   <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "K") %>% 
               select (Time = Time, CK = Conc.)
data_B1_Mu  <- data %>% filter(Study ==7 & Chem == "PFOA" & Species == 'DC' & Matrix == "Mu") %>% 
               select (Time = Time, CM = Conc.)


## Define the fixed parameters
fixpars <- c(
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

## Create a cost function and later used in model optimization  
## Estimate the model residual with experimental data by modcost function (from FME package)
MCost.PFOA<-function (pars){
    
    outdf_B1   <- pred (fixpars=fixpars, pars = pars, BW = 583, dose = 0.002, tdose = 28, atol = 1E-6)
    #outdf_B2   <- pred (fixpars=fixpars, pars = pars, BW = 587, dose = 0.0006/587, tdose = 545, tt = 1)
    
    cost<- modCost  (model = outdf_B1, obs = data_B1_P,  x ="Time", weight = "mean")
    cost<- modCost  (model = outdf_B1, obs = data_B1_L,  x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_B1, obs = data_B1_K,  x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_B1, obs = data_B1_M,  x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_B1, obs = data_B1_Mu, x ="Time", cost = cost, weight = "mean")

    
    return(cost)
}

MCost.PFOA(log(fixpars))

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## Initial and/or optimized parameters
theta_init <- log(c(
    PL           = 1.2,
    PK           = 1.5,
    PM           = 0.0772,
    PU           = 0.2,
    PRest        = 0.05,
    Free         = 0.012,
    KabsC        = 2.12,
    KurineC      = 0.5,
    KbileC       = 0.1,
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
))

MCost.PFOA(theta_init)

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCost.PFOA, parms = theta_init, varscale = 1)
Sen = summary(SnsPlasma)
plot(summary(SnsPlasma))

## Selected sensitive parameters;
theta <- theta_init[abs(Sen$Mean) > 1.5*mean(abs(Sen$Mean))]
theta


##
theta <- log(c(
    #PL           = 1.2,
    #PK           = 1.5,
    #PM           = 0.0772,
    #PU           = 0.2,
    #PRest        = 0.05,
    Free         = 0.012,
    #KabsC        = 2.12,
    #KurineC      = 0.5,
    KbileC       = 0.06,
    KehcC        = 0.002,
    #Kdif         = 0.1,
    KeffluxC     = 0.02,
    KfecesC      = 0.06,
    #RAF_apical   = 0.0005,
    Vmax_apical_invitro = 0.947
    #Km_apical    = 52.3    
    #RAF_baso     = 4.07,
    #Vmax_baso_invitro = 0.4,
    #Km_baso      = 27
    #PMilkM       = 0.9587
    #Kmilking     = 0.5801       
))


# PBPK model fitting 
## Least squares fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
Fit_PFOA <- modFit(f = MCost.PFOA, p = theta, method ="Marq",
              control = nls.lm.control(nprint = 1))

summary(Fit_PFOA)                                  ## Summary of fit 
exp(Fit_PFOA$par)                                  ## Get the arithmetic value out of the log domain

MCost.PFOA(Fit_PFOA$par)


## 
DC_PFOA<-MCost.PFOA(Fit_PFOA$par)


PDat <- cbind.data.frame (OBS = DC_PFOA$residuals$obs,
                          PRE = DC_PFOA$residuals$mod,
                          RES = DC_PFOA$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "DC")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat)
summary(fit)

PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,4), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,4),labels = scales::math_format(10^.x))


p1

## Time vs. concentrations plot for PFOA----------------------------------------
new_pars <- Fit_PFOA$par
DC_outdf_B1   <- pred (fixpars=fixpars, pars = new_pars, BW = 583, dose = 0.002, tdose = 28, endtime = 40*24)

DC_PFOA_pre <- rbind.data.frame(
    
    Pred_B1_P  <- DC_outdf_B1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study    = "B1_P",  Chemical = 'PFOA'),
    Pred_B1_L  <- DC_outdf_B1  %>% select(Time = Time, Conc. = CL)  %>% mutate(Study    = "B1_L",  Chemical = 'PFOA'),
    Pred_B1_K  <- DC_outdf_B1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study    = "B1_K",  Chemical = 'PFOA'),
    Pred_B1_M  <- DC_outdf_B1  %>% select(Time = Time, Conc. = CMilk) %>% mutate(Study  = "B1_M",  Chemical = 'PFOA'),
    Pred_B1_MU <- DC_outdf_B1  %>% select(Time = Time, Conc. = CM) %>% mutate(Study    = "B1_Mu",  Chemical = 'PFOA')

)    

DC_OBS_data <- rbind.data.frame(
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



p2 <-
    ggplot(DC_OBS_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = DC_PFOA_pre, aes(Time/24, Conc.))  +
    #annotation_logticks() + 
    # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(scale_y_continuous(trans  = scales::log2_trans(),
    #                                       breaks = scales::trans_breaks("log2", function(x) 2^x),
    #                                       labels = scales::trans_format("log2", scales::math_format(2^.x))))+
    #scale_x_continuous(labels = scales::math_format(10^.x)) +
    facet_wrap (.~Study, scales = "free")


p2


# PFHxS calibration data set  ---------------------------------------------------
# A7: Kowalczyk et al., 2013; Repeated exposure at dose of 0.0046 mg/kg for 28 days; P, M, Mu, l, k


## Study 7: Kowalczyk et al., 2013
data_C1_P   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "P") %>% 
    select (Time = Time, CP = Conc.)
data_C1_M   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "M") %>% 
    select (Time = Time, CMilk = Conc.)
data_C1_L   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "L") %>% 
    select (Time = Time, CL = Conc.)
data_C1_K   <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "K") %>% 
    select (Time = Time, CK = Conc.)
data_C1_Mu  <- data %>% filter(Study ==7 & Chem == "PFHxS" & Species == 'DC' & Matrix == "Mu") %>% 
    select (Time = Time, CM = Conc.)



## Define the fixed parameters
fixpars <- c(
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

## Create a cost function and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
MCost.PFHxS<-function (pars){
    
    outdf_C1   <- pred (fixpars=fixpars, pars = pars, BW = 583, dose = 0.0046, tdose = 28, atol = 1e-4)

    cost<- modCost  (model = outdf_C1, obs = data_C1_P, x ="Time", weight = "mean")
    cost<- modCost  (model = outdf_C1, obs = data_C1_L, x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_C1, obs = data_C1_K, x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_C1, obs = data_C1_M, x ="Time", cost = cost, weight = "mean")
    cost<- modCost  (model = outdf_C1, obs = data_C1_Mu, x ="Time", cost = cost, weight = "mean")

    
    return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## Initial and/or optimized parameters
theta_init <- log(c(
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
  RAF_apical   = 5,
  Vmax_apical_invitro = 0.947,
  Km_apical    = 52.3,    
  RAF_baso     = 4.07,
  Vmax_baso_invitro = 0.04,
  Km_baso      = 27,
  PMilkM       = 0.09587,
  Kmilking     = 0.5801                    
))

MCost.PFHxS(theta_init)

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCost.PFHxS, parms = theta_init, varscale = 1)
Sen = summary(SnsPlasma)
plot(summary(SnsPlasma))

## Selected sensitive parameters;
theta <- theta_init[abs(Sen$Mean) > 1.5*mean(abs(Sen$Mean))]
theta


##
theta <- log(c(
    #PL         = 0.2,
    #PK         = 0.4,
    PM         = 0.0772,
    #PU         = 0.2,
    PRest      = 0.5,
    Free       = 0.052,
    #KabsC      = 2.12,
    KurineC      = 0.02,
    KbileC       = 0.01,
    #KehcC        = 0.01,
    Kdif         = 0.01,
    KeffluxC   = 0.2,
    #KfecesC    = 0.01,
    #RAF_apical   = 5,
    #Vmax_apical_invitro = 0.947,
    Km_apical    = 52.3,    
    #RAF_baso     = 4.07,
    #Vmax_baso_invitro = 0.04
    #Km_baso      = 27,
    PMilkM       = 0.09587
    #Kmilking     = 0.5801                
))


# PBPK model fitting 
## Least squares fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
Fit_PFHXS<- modFit(f = MCost.PFHxS, p = theta, method ="Marq",
                   control = nls.lm.control(nprint = 1))

summary(Fit_PFHXS)                                  ## Summary of fit 
exp(Fit_PFHXS$par)                                  ## Get the arithmetic value out of the log domain

MCost.PFHxS(Fit_PFHXS$par)


##
DC_PFHXS<-MCost.PFHxS(Fit_PFHXS$par)


PDat <- cbind.data.frame (OBS = DC_PFHXS$residuals$obs,
                          PRE = DC_PFHXS$residuals$mod,
                          RES = DC_PFHXS$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "DC")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat)
summary(fit)

PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,4), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,4),labels = scales::math_format(10^.x))


p1


## Time vs. concentrations plot for PFOA----------------------------------------
new_pars <- Fit_PFHXS$par
DC_outdf_C1   <- pred (fixpars=fixpars, pars = new_pars, BW = 583, dose = 0.0046, tdose = 28, atol = 1e-4)

DC_PFHXS_pre <- rbind.data.frame(
    
    Pred_C1_P  <- DC_outdf_C1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study    = "C1_P",  Chemical = 'PFHXS'),
    Pred_C1_L  <- DC_outdf_C1  %>% select(Time = Time, Conc. = CL)  %>% mutate(Study    = "C1_L",  Chemical = 'PFHXS'),
    Pred_C1_K  <- DC_outdf_C1  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study    = "C1_K",  Chemical = 'PFHXS'),
    Pred_C1_M  <- DC_outdf_C1  %>% select(Time = Time, Conc. = CMilk) %>% mutate(Study  = "C1_M",  Chemical = 'PFHXS'),
    Pred_C1_MU <- DC_outdf_C1  %>% select(Time = Time, Conc. = CM) %>% mutate(Study    = "C1_Mu",  Chemical = 'PFHXS')
    
)    

DC_OBS_data <- rbind.data.frame(
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



p2 <-
    ggplot(DC_OBS_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = DC_PFHXS_pre, aes(Time/24, Conc.))  +
    #annotation_logticks() + 
    # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(scale_y_continuous(trans  = scales::log2_trans(),
    #                                       breaks = scales::trans_breaks("log2", function(x) 2^x),
    #                                       labels = scales::trans_format("log2", scales::math_format(2^.x))))+
    #scale_x_continuous(labels = scales::math_format(10^.x)) +
    facet_wrap (.~Study, scales = "free")


p2



# Save the fitting results
saveRDS(Fit_PFOS,  file = 'Fit_DC_PFOS.rds')
saveRDS(Fit_PFOA,  file = 'Fit_DC_PFOA.rds')
saveRDS(Fit_PFHXS, file = 'Fit_DC_PFHXS.rds')

# Save the fitting results
saveRDS(DC_PFOS,  file = 'DC_PFOS.rds')
saveRDS(DC_PFOA,  file = 'DC_PFOA.rds')
saveRDS(DC_PFHXS, file = 'DC_PFHXS.rds')


