## Translated into mrgsolve based on the code published by Chou and Lin (Environ Int, 2019)
library(mrgsolve)   ## for PBPK modeling
library(ggplot2)    ## for graph demonstration
library(dplyr)      ## for dataframe manipulation
library(FME)        ## for sensitive parameter optimization
library(minpack.lm) ## selecting parameter for optimization, combining with FME
library (optimr)

## Loading the R code
source (file = 'Beef_PBPK.R')
#source (file = 'Cow_PBPK.R')

## Loading the mrgsolve model
mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
#mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)

## Read data
data    <- read.csv(file = "Data.csv")


# PFOS calibration data set  ---------------------------------------------------
# A2: Lupton et al., 2014; single oral exposure at dose of 8 mg/kg; P, L, Mu, K, Sp, Lu
# A3: Lupton et al., 2015; single oral exposure at dose of 0.098 and 9.01 mg/kg; P, L, Mu, K
# A5: Drew, 2019, 2021; feeding with water contaminating 0.00042 mg/kg; P, L, Mu, K 
#-------------------------------------------------------------------------------

## Study 1: Guruge et al., 2008
data_A1_P   <- data %>% filter(Study ==1 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.) 

## Study 2: Lupton et al., 2014
data_A2_P   <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.)
data_A2_Mu  <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu")%>% 
               select (Time = Time, CM = Conc.)
data_A2_L   <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "L") %>% 
               select (Time = Time, CL = Conc.)
data_A2_K   <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "K") %>% 
               select (Time = Time, CK = Conc.)
data_A2_U   <- data %>% filter(Study ==2 & Chem == "PFOS" & Species == 'BC' & Matrix == "U") %>% 
               select (Time = Time, CU = Conc.)

## Study 3: Lupton et al., 2015
data_A3_P_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "P" & Dose == 0.098) %>% 
                 select (Time = Time, CP = Conc.)
data_A3_P_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "P" & Dose == 9.1) %>% 
                 select (Time = Time, CP = Conc.)
data_A3_L_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "L" & Dose == 0.098) %>% 
                 select (Time = Time, CL = Conc.)
data_A3_L_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "L" & Dose == 9.1) %>% 
                 select (Time = Time, CL = Conc.)
data_A3_K_1   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "K" & Dose == 0.098) %>% 
                 select (Time = Time, CK = Conc.)
data_A3_K_2   <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "K" & Dose == 9.1) %>% 
                 select (Time = Time, CK = Conc.)
data_A3_Mu_1  <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu" & Dose == 0.098) %>% 
                 select (Time = Time, CM = Conc.)
data_A3_Mu_2  <- data %>% filter(Study ==3 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu" & Dose == 9.1) %>% 
                 select (Time = Time, CM = Conc.)

## Study 5: Drew, 2019, 2021
data_A5_P   <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.)
data_A5_Mu  <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "Mu")%>% 
               select (Time = Time, CM = Conc.)
data_A5_L   <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "L") %>% 
               select (Time = Time, CL = Conc.)
data_A5_K   <- data %>% filter(Study ==5 & Chem == "PFOS" & Species == 'BC' & Matrix == "K") %>% 
               select (Time = Time, CK = Conc.)



## Define the prediction function
pred <- function(fixpars = (mod_bc %>% param), pars, BW, tdose, dose, tt=1, atol = 1E-8, rtol = atol/2, endtime = 24*365*1) {
    
    ## Get out of log domain
    pars <- lapply(pars, exp)           ## Return a list of exp (Parameters for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    BW          = BW                    ## Body weight based on 3-4 yr Belted Galloway cow (https://beltie.org/PDFs/2017/Belted-galloway-Breeders-Manual-Feb-2017.pdf)
    tinterval   = 24                    ## Time interval; 
    TDOSE       = tdose                 ## Total dosing/Dose times; Repeat drinking water dose for 2 yrs
    DOSE        = dose                  ## Input drinking water dose  
    DOSEoral    = DOSE*BW               ## Amount of dose
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    ex.oral <- ev (ID   = 1,            ## One individual
                   time = 0,            ## Dossed start time 
                   amt  = DOSEoral,     ## Amount of dose 
                   ii   = tinterval,    ## Time interval
                   addl = TDOSE - 1,    ## Additional dosing 
                   cmt  = "AST",        ## The dosing compartment: AST Stomach  
                   replicate = FALSE)   ## No replicate

    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, tt) ## Simulation time post 1-yr exposure
    
    ## Simulation of exposure scenario 
    out <- 
        mod_bc %>%  
        param (fixpars) %>%
        param (pars) %>%                 
        update(atol = atol, rtol = rtol, maxsteps = 50000) %>%  
        mrgsim_d (data = ex.oral, tgrid = tsamp) 
    

    outdf <- out %>% as.data.frame
    
    outdf = cbind.data.frame (Time    = outdf$time, 
                              CP      = (outdf$Plasma)*1000, 
                              CL      = (outdf$Liver)*1000,
                              CK      = (outdf$Kidney)*1000, 
                              CM      = (outdf$Muscle)*1000,
                              CU      = outdf$Urine
                              )

    return (outdf) 
}

fixpars <- c(
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

## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
MCost_PFOS<-function (pars, w = 'none'){
    
    outdf_A1   <- pred (fixpars=fixpars, pars = pars, BW = 329, dose = 0.00000798, tdose =30*28, tt = 960)
    outdf_A2   <- pred (fixpars=fixpars, pars = pars, BW = 329, dose = 8,       tdose = 1, tt = 1)
    outdf_A3_1 <- pred (fixpars=fixpars, pars = pars, BW = 343, dose = 0.098,   tdose = 1, tt = 1)
    outdf_A3_2 <- pred (fixpars=fixpars, pars = pars, BW = 343, dose = 9.09,    tdose = 1, tt = 1)
    outdf_A5   <- pred (fixpars=fixpars, pars = pars, BW = 364, dose = 0.00042, tdose = 365*1.5, tt =120)
    
    cost<- modCost  (model = outdf_A1,   obs = data_A1_P,   x ="Time", weight = w)
    cost<- modCost  (model = outdf_A2,   obs = data_A2_P,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A2,   obs = data_A2_Mu,  x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A2,   obs = data_A2_L,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A2,   obs = data_A2_U,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_1, obs = data_A3_P_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_1, obs = data_A3_L_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_1, obs = data_A3_K_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_1, obs = data_A3_Mu_1,x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_2, obs = data_A3_P_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_2, obs = data_A3_L_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_2, obs = data_A3_K_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A3_2, obs = data_A3_Mu_2,x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A5,   obs = data_A5_P,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A5,   obs = data_A5_Mu,  x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A5,   obs = data_A5_L,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = outdf_A5,   obs = data_A5_K,   x ="Time", cost = cost, weight = w)

    return(cost)
}


## Local sensitivity analysis
## Choose the senstiive parameters in the model
## Optimized parameters
theta.init <- log(c(
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
))

MCost_PFOS(theta.init, w = 'mean')

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
# SnsPlasma <- sensFun(func = MCost_PFOS, parms = theta.init, varscale = 1)
# Sen = summary(SnsPlasma)
# plot(summary(SnsPlasma))
# 

## Selected sensitive parameters;
# theta <- theta.init[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
# theta

##
theta <- log(c(
    #PL           = 1.03,
    #PK           = 0.48,
    PM           = 0.08,
    PRest        = 0.15,
    Free         = 0.055,
    #KabsC        = 2.12,
    KurineC      = 0.01,
    KbileC       = 0.01,
    KehcC        = 0.01,
    #Kdif         = 0.01,
    #KeffluxC     = 0.1,
    #KfecesC      = 0.01,
    RAF_apical   = 0.001,
    #Vmax_apical_invitro = 51803,
    Km_apical    = 64.400
    #RAF_baso     = 0.005
    #Vmax_baso_invitro = 60.99,
    #Km_baso      = 34.500
))



## PBPK model fitting 
## Least squares fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
Fit_PFOS <- modFit(f = MCost_PFOS, w = 'mean', p = theta, method = "Marq",
               control = nls.lm.control(nprint = 1))



summary(Fit_PFOS)  ## Summary of fit 
exp(Fit_PFOS$par)                                  ## Get the arithmetic value out of the log domain

MCost_PFOS(Fit_PFOS$par, w ='mean')


## Global evaluation of goodness of fit------------------------------------------
BC_PFOS<-MCost_PFOS(Fit_PFOS$par)


PDat_PFOS <- cbind.data.frame (OBS = BC_PFOS$residuals$obs,
                          PRE = BC_PFOS$residuals$mod,
                          RES = BC_PFOS$residuals$res)

PDat_PFOS <- PDat_PFOS %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Chemical = "PFOS")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat_PFOS)
summary(fit)

PlotDat_PFOS <- PDat_PFOS %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat_PFOS, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Chemical)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,5), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-2,5),labels = scales::math_format(10^.x))


p1

## Time vs. concentrations plot ------------------------------------------------

new_pars <- Fit_PFOS$par
BC_outdf_A1   <- pred (fixpars=fixpars, pars = new_pars, BW = 329, dose = 0.00000798, tdose =30*28, tt = 960)
BC_outdf_A2   <- pred (fixpars=fixpars, pars = new_pars, BW = 329, dose = 8,       tdose = 1, tt = 1)
BC_outdf_A3_1 <- pred (fixpars=fixpars, pars = new_pars, BW = 343, dose = 0.098,   tdose = 1, tt = 1)
BC_outdf_A3_2 <- pred (fixpars=fixpars, pars = new_pars, BW = 343, dose = 9.09,    tdose = 1, tt = 1)
BC_outdf_A5   <- pred (fixpars=fixpars, pars = new_pars, BW = 364, dose = 0.00042, tdose = 365*1.5, tt =120)


BC_PFOS_pre <- rbind.data.frame(
    
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
    Pred_A3_L_1   <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CL) %>% mutate(Study = "A3_L_1",  Chemical = 'PFOS'),
    Pred_A3_L_2   <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CL) %>% mutate(Study = "A3_L_2",  Chemical = 'PFOS'),
    Pred_A3_K_1   <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A3_K_1",  Chemical = 'PFOS'),
    Pred_A3_K_2   <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A3_K_2",  Chemical = 'PFOS'),
    Pred_A3_Mu_1  <- BC_outdf_A3_1 %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A3_Mu_1", Chemical = 'PFOS'),
    Pred_A3_Mu_2  <- BC_outdf_A3_2 %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A3_Mu_2", Chemical = 'PFOS'),
    
    ## Study 5: Drew, 2019, 2021
    Pred_A5_P     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CP) %>% mutate(Study = "A5_P", Chemical = 'PFOS'),
    Pred_A5_Mu    <- BC_outdf_A5  %>% select (Time = Time, Conc. = CM) %>% mutate(Study = "A5_Mu", Chemical = 'PFOS'),
    Pred_A5_L     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CL) %>% mutate(Study = "A5_L", Chemical = 'PFOS'),
    Pred_A5_K     <- BC_outdf_A5  %>% select (Time = Time, Conc. = CK) %>% mutate(Study = "A5_K", Chemical = 'PFOS')
)


BC_OBS_data <- rbind.data.frame(
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
    
    
    

p2 <-
    ggplot(BC_OBS_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = BC_PFOS_pre, aes(Time/24, Conc.))  +
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
# A4: Lupton et al., 2012; single oral exposure at dose of 8 mg/kg; P, U, F
#-------------------------------------------------------------------------------

data_B1_P   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "P") %>% 
               select (Time = Time, CP = Conc.)
data_B1_U   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "U") %>% 
               select (Time = Time, CU = Conc.)



## Define the fixed parameters
fixpars <- c(
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




## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
MCost_PFOA<-function (pars){
    
  
    outdf_B1   <- pred (pars=pars, BW = 329, dose = 1, tdose = 1, tt = 1, atol = 1E-8, rtol = 1E-4)

    cost<- modCost  (model = outdf_B1, obs = data_B1_P, x = "Time",weight = "mean")
    cost<- modCost  (model = outdf_B1, obs = data_B1_U, cost = cost, x ="Time", weight = "mean")
    
    return(cost)
}

##
theta.init <- log(c(
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
))


MCost_PFOA(pars=theta.init)



## Sensitivity function (FME) 
## Check the sensitive parameters in the model
SnsPlasma <- sensFun(func = MCost_PFOA, parms = theta.init, varscale = 1)
Sen = summary(SnsPlasma)
plot(summary(SnsPlasma))


## Selected sensitive parameters;
theta <- theta.init[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
theta


##
theta <- log(c(
    #PL           = 1.03,
    #PK           = 0.48,
    #PM           = 0.08,
    PRest        = 0.1526,
    Free         = 0.055,
    #KabsC        = 2.12,
    KurineC      = 0.1,
    #KbileC       = 0.01
    #KehcC        = 0.001,
    #Kdif         = 0.001,
    #KeffluxC     = 0.1,
    #KfecesC      = 0.0137,
    #RAF_apical   = 1e-3,
    Vmax_apical_invitro = 0.9
    #Km_apical    = 64.400,
    #RAF_baso     = 0.005,
    #Vmax_baso_invitro = 6.9
    #Km_baso      = 34.500
))


## PBPK model fitting 
## Least squares fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
Fit_PFOA <- modFit(f = MCost_PFOA, p = theta, method ="Marq",
              control = nls.lm.control(nprint = 1))

summary(Fit_PFOA)                                  ## Summary of fit 
exp(Fit_PFOA$par)                                  ## Get the arithmetic value out of the log domain

MCost_PFOA(Fit_PFOA$par)



##
BC_PFOA<-MCost_PFOA(Fit_PFOA$par)


PDat_PFOA <- cbind.data.frame (OBS = BC_PFOA$residuals$obs,
                          PRE = BC_PFOA$residuals$mod,
                          RES = BC_PFOA$residuals$res)

PDat_PFOA <- PDat_PFOA %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Chmeical = "PFOA")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat_PFOA)
summary(fit)

PlotDat_PFOA <- PDat_PFOA %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat_PFOA, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Chmeical)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(1,4), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(1,4), labels = scales::math_format(10^.x))


p1

### Time vs. concentration plot ------------------------------------------------

new_pars <- Fit_PFOA$par
BC_outdf_B1   <- pred (fixpars=fixpars, pars=new_pars, BW = 329, dose = 1, tdose = 1, tt = 1, atol = 1E-6, rtol = 1E-3, endtime = 24*10)


BC_PFOA_pre <- rbind.data.frame(
    
    ## Study 1: Guruge et al., 2008
    Pred_B1_P     <- BC_outdf_B1  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "B1_P",  Chemical = 'PFOA'),
    Pred_B1_U     <- BC_outdf_B1  %>% select(Time = Time, Conc. = CU)  %>% mutate(Study = "B1_U",  Chemical = 'PFOA')
    
)


BC_PFOA_data <- rbind.data.frame(
    data_B1_P   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "P") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_P", Chemical = 'PFOA'),
    data_B1_U   <- data %>% filter(Study ==4 & Chem == "PFOA" & Species == 'BC' & Matrix == "U") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "B1_U", Chemical = 'PFOA')
)




p2 <-
    ggplot(BC_PFOA_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = BC_PFOA_pre, aes(Time/24, Conc.))  +
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
# A5: Drew, 2019, 2021; 
#-------------------------------------------------------------------------------
data_C1_P   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "P" & Time > 15120) %>% 
               select (Time = Time, CP = Conc.)
data_C2_L   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "L") %>% 
               select (Time = Time, CL = Conc.)
data_C3_K   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "K") %>% 
               select (Time = Time, CK = Conc.)
data_C4_Mu  <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "Mu") %>% 
               select (Time = Time, CM = Conc.)


## Define the fixed parameters
fixpars <- c(
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


## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
MCost_PHXS<-function (pars){
    
    
    outdf_C   <- pred (fixpars=fixpars, pars, BW = 364, dose = 0.00084, tdose = 365*2, tt =120)
    
    cost<- modCost  (model = outdf_C, obs = data_C1_P, x = "Time",weight = "mean")
    cost<- modCost  (model = outdf_C, obs = data_C2_L,  cost = cost, x = "Time",  weight = "mean")
    cost<- modCost  (model = outdf_C, obs = data_C3_K,  cost = cost, x = "Time", weight = "mean")
    cost<- modCost  (model = outdf_C, obs = data_C4_Mu, cost = cost, x = "Time",  weight = "mean")
    
    return(cost)
}

##
theta.init <- log(c(
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
))


MCost_PHXS(theta.init)

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
# SnsPlasma <- sensFun(func = MCost_PHXS, parms = theta.init, varscale = 1)
# Sen = summary(SnsPlasma)
# plot(summary(SnsPlasma))


## Selected sensitive parameters;
theta <- theta.init[abs(Sen$Min) > 1.2*mean(abs(Sen$Min))]
theta

##
theta <- log(c(
    #PL           = 0.21,
    #PK           = 0.27,
    #PM           = 0.062,
    PRest        = 0.0715,
    Free         = 0.02,
    #KabsC        = 2.21,
    #KurineC      = 0.01,
    #KbileC       = 0.1,
    #KehcC        = 0.01,
    #Kdif         = 0.001,
    #KeffluxC     = 0.1,
    #KfecesC      = 0.0137,
    RAF_apical   = 1e-4,
    Vmax_apical_invitro = 51803
    #Km_apical    = 64.400    
    #RAF_baso     = 1e-3,
    #Vmax_baso_invitro = 60.99
    #Km_baso      = 34.500
))

## training gestational parameters
Fit_PFHXS <- modFit(f = MCost_PHXS, p = theta, method ="Marq",
              control = nls.lm.control(nprint = 1))

summary(Fit_PFHXS)                                  ## Summary of fit 
exp(Fit_PFHXS$par)                                  ## Get the arithmetic value out of the log domain

MCost_PHXS(Fit_PFHXS$par)



##
BC_PFHXS<-MCost_PHXS(Fit_PFHXS$par)


PDat_PFHXS <- cbind.data.frame (OBS = BC_PFHXS$residuals$obs,
                          PRE = BC_PFHXS$residuals$mod,
                          RES = BC_PFHXS$residuals$res)

PDat_PFHXS <- PDat_PFHXS %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Chmeical = "PFHXS")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat_PFHXS)
summary(fit)

## PLOT TO CHECK THE GOODNESS OF FIT
PlotDat_PFHXS <- PDat_PFHXS %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat_PFHXS, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Chmeical)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x))


p1

### Time vs. concentration plot ------------------------------------------------

new_pars <- Fit_PFHXS$par
BC_outdf_C   <- pred (fixpars=fixpars, new_pars, BW = 364, dose = 0.00084, tdose = 365*2, tt =120)


BC_PFHxS_pre <- rbind.data.frame(
    
    Pred_C1_P     <- BC_outdf_C  %>% select(Time = Time, Conc. = CP)  %>% mutate(Study = "C1_P",  Chemical = 'PFHxS'),
    Pred_C1_L     <- BC_outdf_C  %>% select(Time = Time, Conc. = CL)  %>% mutate(Study = "C1_L",  Chemical = 'PFHxS'),
    Pred_C1_K     <- BC_outdf_C  %>% select(Time = Time, Conc. = CK)  %>% mutate(Study = "C1_K",  Chemical = 'PFHxS'),
    Pred_C1_M     <- BC_outdf_C  %>% select(Time = Time, Conc. = CM)  %>% mutate(Study = "C1_M",  Chemical = 'PFHxS')
    
)


BC_PFHxS_data <- rbind.data.frame(
    
    data_C1_P   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "P" & Time > 15120) %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_P", Chemical = 'PFHxS'),
    data_C1_L   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "L") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_L", Chemical = 'PFHxS'),
    data_C1_K   <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "K") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_K", Chemical = 'PFHxS'),
    data_C1_Mu  <- data %>% filter(Study ==5 & Chem == "PFHxS" & Species == 'BC' & Matrix == "Mu") %>% 
        select (Time = Time, Conc. = Conc.) %>% mutate(Study = "C1_M", Chemical = 'PFHxS')
    
)




p2 <-
    ggplot(BC_PFHxS_data, aes(Time/24, Conc.)) +
    geom_point  ( colour = "steelblue4", size = 1)  +
    geom_line  (data = BC_PFHxS_pre, aes(Time/24, Conc.))  +
    #annotation_logticks() + 
    # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(scale_y_continuous(trans  = scales::log2_trans(),
    #                                       breaks = scales::trans_breaks("log2", function(x) 2^x),
    #                                       labels = scales::trans_format("log2", scales::math_format(2^.x))))+
    scale_x_continuous(limit = c(730, 900)) +
    facet_wrap (.~Study, scales = "free")


p2


# Save the fitting results------------------------------------------------------
saveRDS(Fit_PFOS,  file = 'Fit_BC_PFOS.rds')
saveRDS(Fit_PFOA,  file = 'Fit_BC_PFOA.rds')
saveRDS(Fit_PFHXS, file = 'Fit_BC_PFHXS.rds')

# Save the fitting results
saveRDS(BC_PFOS,  file = 'BC_PFOS.rds')
saveRDS(BC_PFOA,  file = 'BC_PFOA.rds')
saveRDS(BC_PFHXS, file = 'BC_PFHXS.rds')



