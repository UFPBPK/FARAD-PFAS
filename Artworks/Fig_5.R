
## MODEL EVALUATION (PFOA and PFOS for beef cattle)-----------------------------
# - Author        : Wei-Chun Chou
# - Supervisor    : Zhoumeng Lin
# - Date          : Started in February 2022, finalized in June 2023
# - Abbreviation  : BC, Beef cattle; DC, Dairy cow  
#------------------------------------------------------------------------------

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
source (file = 'Function.R')  # Loading the parameters


## Loading the mrgsolve model
mod_bc   <- mcode_cache ("PFOS_Beef_PBPK.code", PFOS_Beef_PBPK.code)
mod_dc   <- mcode_cache ("PFOS_Cow_PBPK.code", PFOS_Cow_PBPK.code)


# Define the mrg prediction function
### Define the prediction function
pred <- function (pars, Species, Dose, idata, route = 'oral', tinterval = 24, Dtimes = 1, tt = 24, N = 10, endtime = 24*30) {
    
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
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + endtime, tt) 
    
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
                       SD_Milk   = sd(Milk)
                       
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

# - Water samples in Xinjiang area (ref.6: Table 1): PFOA (0.75 ng/L), PFOS (0.34 ng/L), PFHXS (0.25) 
# - Soil samples: PFOA: 354 (ng/kg); PFOS: 193 (ng/kg); PFHXS: 27.9 (ng/kg)

# Reference list-------------
# Ref.1: Maine CDC (https://www.maine.gov/dep/spills/topics/pfas/Agronomic-Pathway-Soil-Screening-Levels-Soil-Fodder-Cows-Milk-09.16.20.pdf)
# Ref.2: EPA PRG (https://epa-prgs.ornl.gov/radionuclides/)
# Ref.3: Beef Cattle Water Requirements and Source Management (https://reurl.cc/Zr2DRl)
# Ref.4: Wang et al. (2017); doi: 10.3390/ijerph14090970
# Ref.5: Xing et al. (2016); doi: 10.3390/ijerph13101037
# Ref.6: Li et al. (2019): doi: https://doi.org/10.1016/j.envint.2018.11.036
# Ref.7: Li et al. (2020): doi: https://doi.org/10.1016/j.envint.2019.105419
#//////////////////////////////////////////////////////////////////////////////

## Calculation of idata
idata_BC_PFOS <- MCsim (pars_BC_PFOS, pars_SD_BC)
idata_BC_PFOA <- MCsim (pars_BC_PFOA, pars_SD_BC)

idata_DC_PFOS <- MCsim (pars_DC_PFOS, pars_SD_DC)
idata_DC_PFOA <- MCsim (pars_DC_PFOA, pars_SD_DC)


## Cattle
Sim_PFOA_BC <-pred (pars = pars_BC_PFOA, Species = 'BC', Dose = (5.19E-6), Dtimes = 365*2, tt= 24,  idata = idata_BC_PFOA)
Sim_PFOS_BC <-pred (pars = pars_BC_PFOS, Species = 'BC', Dose = (6.26E-7), Dtimes = 365*2, tt= 24,  idata = idata_BC_PFOS)

Sim_PFOA_DC <-pred (pars = pars_DC_PFOA, Species = 'DC', Dose = (5.88E-6), Dtimes = 365*3, tt= 24,  idata = idata_DC_PFOA)
Sim_PFOS_DC <-pred (pars = pars_DC_PFOS, Species = 'DC', Dose = (7.63e-7), Dtimes = 365*3, tt= 24,  idata = idata_DC_PFOS)

data_pfos_bc <- Sim_PFOS_BC %>% filter (row_number() == 365*2)
data_pfoa_bc <- Sim_PFOA_BC %>% filter (row_number() == 365*2)
data_pfos_dc <- Sim_PFOS_DC %>% filter (row_number() == 365*3)
data_pfoa_dc <- Sim_PFOA_DC %>% filter (row_number() == 365*3)





Pdata_bc <- cbind.data.frame(  Mean  = c(data_pfos_bc$M_Muscle*1000, data_pfos_bc$M_Liver*1000,
                                           data_pfoa_bc$M_Muscle*1000, data_pfoa_bc$M_Liver*1000,
                                           data_pfoa_dc$M_Milk*1000, data_pfos_dc$M_Milk*1000, 0.008, 0.617, 0.011, 0.264, 16.2e-3, 24.5e-3),
                            
                               SD      = c(data_pfos_bc$SD_Muscle*1000, data_pfos_bc$SD_Liver*1000,
                                           data_pfoa_bc$SD_Muscle*1000, data_pfoa_bc$SD_Liver*1000,
                                           data_pfoa_dc$SD_Milk*1000, data_pfoa_dc$SD_Milk*1000, 0.0042, 0.27, 0.00098, 0.1246, 0, 0),
                               Type     = c(rep("pred",6), rep("obs",6)),
                               Target   = c(rep(c("PFOS-Muscle", "PFOS-Liver", "PFOA-Muscle", "PFOA-Liver", "PFOA-Milk", "PFOS-Milk"),2)),
                               Chemical = rep(c("PFOS", "PFOS", "PFOA", "PFOA", "PFOA", "PFOS"),2),
                               Species  = c(rep('BC',4), rep('DC',2),rep('BC',4), rep('DC',2))
               
)

                 
Pdata_ratio <- cbind.data.frame(OP_ratio = c(0.008/0.045, 0.617/0.97, 0.01/0.01, 0.26/0.08, 0.016/0.01, 0.0245/0.022),
                                 Target   = c("PFOS-Muscle", "PFOS-Liver", "PFOA-Muscle", "PFOA-Liver", "PFOA-Milk", "PFOS-Milk"))            
       
# Change the colors manually
p1<-ggplot(Pdata_bc %>%filter(Target!="PFOS-Liver" & Target!="PFOA-Liver")) + 
    geom_bar(aes(Target , Mean, fill=Type),stat="identity", color="black", position=position_dodge())+
    geom_errorbar(aes(x=Target , y=Mean, ymin=Mean, ymax=Mean+SD, group=Type), width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(limits = c(0,0.15))+
    #geom_point(data=Pdata_ratio, aes(factor(Target), OP_ratio))+
    #scale_y_continuous(sec.axis = sec_axis(~.*1.2))+
    labs (x = "", y = "") +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, linewidth=2),
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
ggsave("Fig 5-1.tiff",
       plot = p1,
       width = 30, height = 20, units = "cm", dpi=320)





p2<- ggplot(Pdata_bc %>%filter(Target %in% c("PFOS-Liver","PFOA-Liver"))) + 
    geom_bar(aes(Target , Mean, fill=Type),stat="identity", color="black", position=position_dodge())+
    geom_errorbar(aes(x=Target , y=Mean, ymin=Mean, ymax=Mean+SD, group=Type), width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(limits = c(0,2))+
    labs (x = "", y = "") +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, linewidth=2),
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
ggsave("Fig 5-2.tiff",
       plot = p2,
       width = 20, height = 10, units = "cm", dpi=320)



# CALCULATE R2
a <- cbind.data.frame(
A <- Pdata_bc %>%filter(Type=='pred') %>% select(Mean),
B  <- Pdata_bc %>%filter(Type=='obs') %>% select(Mean)
)

colnames(a) <- c('pred','obs')

fit_r <- lm(as.numeric(obs) ~ as.numeric(pred), data=a)
summary(fit_r)







