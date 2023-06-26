#------------------------------------------------------------------------------
# This code is used to run the global goodness of fit plot
#------------------------------------------------------------------------------

## Loading R packages
library(ggplot2)
library(dplyr)
library(ggprism) 
library(patchwork)
library(DescTools)
library(ggExtra)
library(gridExtra)

## Loading RDS files; These RDS files produced from model calibration
BC_PFOS  <- readRDS(file = 'BC_PFOS.rds')
BC_PFOA  <- readRDS(file = 'BC_PFOA.rds')
BC_PFHXS <- readRDS(file = 'BC_PFHXS.rds')

DC_PFOS  <- readRDS(file = 'DC_PFOS.rds')
DC_PFOA  <- readRDS(file = 'DC_PFOA.rds')
DC_PFHXS <- readRDS(file = 'DC_PFHXS.rds')

## Make the table for plot the figure for 'PFOS'
## plot data of cattle 
PDat_BC_PFOS <- cbind.data.frame (OBS = BC_PFOS$residuals$obs, 
                                  PRE = BC_PFOS$residuals$mod,
                                  RES = BC_PFOS$residuals$res,
                                  Species = "BC",
                                  Chemical = "PFOS")

PDat_BC_PFOA <- cbind.data.frame (OBS = BC_PFOA$residuals$obs, 
                                  PRE = BC_PFOA$residuals$mod,
                                  RES = BC_PFOA$residuals$res,
                                  Species = "BC",
                                  Chemical = "PFOA")

PDat_BC_PFHXS <- cbind.data.frame (OBS = BC_PFHXS$residuals$obs, 
                                  PRE = BC_PFHXS$residuals$mod,
                                  RES = BC_PFHXS$residuals$res,
                                  Species = "BC",
                                  Chemical = "PFHXS")

## plot data of cow 
PDat_DC_PFOS <- cbind.data.frame (OBS = DC_PFOS$residuals$obs, 
                                  PRE = DC_PFOS$residuals$mod,
                                  RES = DC_PFOS$residuals$res,
                                  Species = "DC",
                                  Chemical = "PFOS")

PDat_DC_PFOA <- cbind.data.frame (OBS = DC_PFOA$residuals$obs, 
                                  PRE = DC_PFOA$residuals$mod,
                                  RES = DC_PFOA$residuals$res,
                                  Species = "DC",
                                  Chemical = "PFOA")

PDat_DC_PFHXS <- cbind.data.frame (OBS = DC_PFHXS$residuals$obs, 
                                   PRE = DC_PFHXS$residuals$mod,
                                   RES = DC_PFHXS$residuals$res,
                                   Species = "DC",
                                   Chemical = "PFHXS")

## combine the data of cattle and swine
PDat_BC <- rbind(PDat_BC_PFOS, PDat_BC_PFOA,PDat_BC_PFHXS)
PDat_DC <- rbind(PDat_DC_PFOS, PDat_DC_PFOA,PDat_DC_PFHXS)

PDat_BC <- PDat_BC %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))
PDat_DC <- PDat_DC %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))


## Estimate the R-squared using linear regression
fit_BC <- lm(Log.PRE ~ Log.OBS, data = PDat_BC)
summary(fit_BC)

fit_DC <- lm(Log.PRE ~ Log.OBS, data = PDat_DC)
summary(fit_DC)


## Add the observed-to-prediction ratios column
PlotDat_BC <- PDat_BC %>% mutate(prediction = predict(fit_BC), OPR = PRE/OBS)
PlotDat_DC <- PDat_DC %>% mutate(prediction = predict(fit_DC), OPR = PRE/OBS)

## Plot theme
ptheme<-theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none') 


## GGplot

## Plot for model calibration
PlotDat <- rbind(PlotDat_BC,PlotDat_DC)
p1 <- 
    ggplot(PlotDat_BC, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Chemical),  shape = factor(Species)), size = 4)+
    scale_colour_brewer(palette = "Dark2") +
    #scale_alpha_manual(values = c(.7, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,5), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,5),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")


p1


p2 <- 
    ggplot(PlotDat_DC, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Chemical),  shape = factor(Species)), size = 4)+
    scale_colour_brewer(palette = "Dark2") +
    #scale_alpha_manual(values = c(.7, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,5), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,5),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")


p2


p1+p2


#//////////////////////////////////////////////////////////////////////////////


p2 <-
    ggplot(PlotDat, aes(Log.PRE, log(OPR,10))) +
    geom_hline(yintercept = log10(2),linetype = 3,color   = "red", size =1) +
    geom_hline(yintercept = log10(0.5),linetype = 3,color   = "red", size =1) +
    #geom_ref_line(colour = "white", h = 0, size =0.5) +
    geom_point(aes(colour = factor(Chemical),  shape = factor(Species)), size = 4) +
    geom_smooth(se = FALSE) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,5),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")


ggMarginal(p2, type = "histogram", margins = "y",  
           yparams = list(binwidth = 0.1, fill = "#FFA488"))


p1+ggMarginal(p2, type = "histogram", margins = "y",  
           yparams = list(binwidth = 0.1, fill = "#FFA488"))











