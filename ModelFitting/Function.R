## Loading R-packages
library(ggplot2)
library(extrafont)
library(ggh4x)
library(EnvStats)
library(dplyr)

## Input the rds file
Fit_BC_PFOS  = readRDS(file = 'Fit_BC_PFOS.rds')
Fit_BC_PFOA  = readRDS(file = 'Fit_BC_PFOA.rds')
Fit_BC_PFHXS = readRDS(file = 'Fit_BC_PFHXS.rds')
Fit_DC_PFOS  = readRDS(file = 'Fit_DC_PFOS.rds')
Fit_DC_PFOA  = readRDS(file = 'Fit_DC_PFOA.rds')
Fit_DC_PFHXS = readRDS(file = 'Fit_DC_PFHXS.rds')


## Define the values of Beef cattle (BC) physiological parameters
## Source: Table 2 and Table 21 (Lin et al., 2020)
Phy_pars_BC <- c(
    BW          = 419,  ## Male beef cattle (Lin et al., 2020)
    QCC         = 5.45, ## Table 13 (Lin et al., 2020)   
    QLCa        = 0.44,	              
    QKCa        = 0.11,	              
    QMCa        = 0.28,	                      
    QRestCa     = (1-0.44-0.11-0.28),
    VLCa        = 0.012,                                     
    VKCa        = 0.002,                                  
    VMCa        = 0.36,                                  
    VbloodCa    = 0.039,
    VRestCa     = (1-0.012-0.002-0.36-0.039)
)

## Define the SD of  Beef cattle (BC) physiological parameters
## Source: Table 2 and Table 21 (Lin et al., 2020)
pars_SD_BC <- c(
    BW          = 46.9, ## Male beef cattle (Lin et al., 2020)
    QCC         = 1.47, ## Table 13 (Lin et al., 2020)    
    QLCa        = 0.25,	              
    QKCa        = 0.08,	              
    QMCa        = 0.09,	                      
    QRestCa     = (1-0.44-0.11-0.28)*0.3, ## assumed CV = 30%
    VLCa        = 0.0018,                                     
    VKCa        = 0.0005,                                  
    VMCa        = 0.12,                                  
    VbloodCa    = 0.0068,
    VRestCa     = (1-0.012-0.002-0.36-0.039)*0.3 ## assumed CV = 30%
)

## Define the values of diary cow (DC) physiological parameters
## Source: Table 4 and Table 23 from (Lin et al., 2020) and table from (Li et al. 2018) 
Phy_pars_DC <- c(
    BW          = 583, ## values collected from Kowalczyk et al. (2013)
    QCC         = 5.45,    
    QLCa        = 0.58,	              
    QKCa        = 0.07,	              
    QMCa        = 0.18,
    QUCa        = 0.13,
    QRestCa     = (1-0.58-0.07-0.18-0.13),
    VLCa        = 0.0125,                                     
    VKCa        = 0.0022,                                  
    VMCa        = 0.36,   # (Li et al. 2018)    
    VUCa        = 0.0318,
    VbloodCa    = 0.0433,
    VRestCa     = (1-0.0135-0.0022-0.27-0.0314-0.0433)
)

# Define the SD of diary cow (DC)  physiological parameters
## Source: Table 4 and Table 23 from (Lin et al., 2020) and table from (Li et al. 2018) 
pars_SD_DC <- c(
    BW          = 31,   ## values collected from Kowalczyk et al. (2013)
    QCC         = 1.47, ## Table 13 (Lin et al., 2020)    
    QLCa        = 0.19,	              
    QKCa        = 0.0021,	              
    QMCa        = 0.054,
    QUCa        = 0.006,
    QRestCa     = (1-0.58-0.07-0.18-0.13)*0.3,
    VLCa        = 0.0028,                                     
    VKCa        = 0.0004,                                  
    VMCa        = 0.081,   ## assumed CV = 30%                              
    VUCa        = 0.0071,
    VbloodCa    = 0.0088,
    VRestCa     = (1-0.0135-0.0022-0.27-0.0314-0.0433)*0.3 ## assumed CV = 30%
)

## Chemical parameters for beef cattle
Che_BC_PFOS <- c(
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

Che_BC_PFOA <- c(
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


Che_BC_PFHXS <- c(
    PL           = 0.21,
    PK           = 0.27,
    PM           = 0.062,
    PRest        = 0.0715,
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

##------------------------------------------------------------------------------
## Chemical parameters for Dairy cow
Che_DC_PFOS <- c(
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

Che_DC_PFOA <- c(
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
    Vmax_apical_invitro = 0.947,
    Km_apical    = 52.3,    
    RAF_baso     = 4.07,
    Vmax_baso_invitro = 0.4,
    Km_baso      = 27,
    PMilkM       = 0.9587,
    Kmilking     = 0.5801       
)


Che_DC_PFHXS <- c(
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

##////////////////////////////////////////////////////////////////////////////
## Define the chemical parameters
## Extract the best optimized parameters and instead the initial parameters with them
Che_BC_PFOS_fit  <- exp(Fit_BC_PFOS$par)
Che_BC_PFOA_fit  <- exp(Fit_BC_PFOA$par)
Che_BC_PFHXS_fit <- exp(Fit_BC_PFHXS$par)


Che_BC_PFOS  [names(Che_BC_PFOS_fit)]   <- as.numeric(Che_BC_PFOS_fit) 
Che_BC_PFOA  [names(Che_BC_PFOA_fit)]   <- as.numeric(Che_BC_PFOA_fit) 
Che_BC_PFHXS [names(Che_BC_PFHXS_fit)]  <- as.numeric(Che_BC_PFHXS_fit) 


## Define dairy cow chemical parameters
Che_DC_PFOS_fit  <- exp(Fit_DC_PFOS$par)
Che_DC_PFOA_fit  <- exp(Fit_DC_PFOA$par)
Che_DC_PFHXS_fit <- exp(Fit_DC_PFHXS$par)


Che_DC_PFOS [names(Che_DC_PFOS_fit)]    <- as.numeric(Che_DC_PFOS_fit) 
Che_DC_PFOA [names(Che_DC_PFOA_fit)]    <- as.numeric(Che_DC_PFOA_fit) 
Che_DC_PFHXS [names(Che_DC_PFHXS_fit)]  <- as.numeric(Che_DC_PFHXS_fit) 


##
pars_BC_PFOS  <- c(Phy_pars_BC, Che_BC_PFOS)
pars_BC_PFOA  <- c(Phy_pars_BC, Che_BC_PFOA)
pars_BC_PFHXS <- c(Phy_pars_BC, Che_BC_PFHXS)

pars_DC_PFOS  <- c(Phy_pars_DC, Che_DC_PFOS)
pars_DC_PFOA  <- c(Phy_pars_DC, Che_DC_PFOA)
pars_DC_PFHXS <- c(Phy_pars_DC, Che_DC_PFHXS)


## Define the Monte Carlo simulation function to output "idata"
MCsim <- function(pars, sd_pars, N=1000) {
    
    
    ## Define the vector
    CV   <-vector()
    M    <- vector()
    Dist <- vector()
    
    ## Make a loop
    for (i in 1:length(pars)){
        
        CV[i]  = switch(substr(names(pars)[i],1,1),
                        "F"   = 0.1,   ## CV of percentage of plasma protein: 0.1 
                        "B"   = 0.2,   ## CV of BW 
                        "H"   = 0.1,   ## CV of htc: 0.1
                        "Q"   = 0.3,#as.numeric(sd_pars[i])/as.numeric(pars[i]),   ## CV of blood flow: 0.3
                        "V"   = 0.3,#if_else(i<10, as.numeric(sd_pars[i])/as.numeric(pars[i]), 0.3),   ## CV of tissue volume: 0.3
                        "K"   = 0.2,   ## CV of kinetic constant: 0.3
                        "P"   = 0.2,   ## CV of partition coefficient: 0.2 
                        0.1)
        
        Dist[i] = switch(substr(names(pars)[i],1,1),
                         "F"   = 2, ## 1 indicate the normal distribution
                         "B"   = 1, ## 2 indicate the log-normal distribution
                         "H"   = 1,
                         "Q"   = 2,
                         "V"   = 1,
                         "K"   = 2,
                         "P"   = 2,
                         1)
        
        M[i] = as.numeric(pars[i])
        
    }
    
    Dat <- cbind.data.frame(Mean = M, CV = CV, Dist = Dist) %>% 
        mutate(SD    = M*CV, 
               logM  = log(M/sqrt(1+(SD/M)^2)),
               logSD = sqrt(log(1+(SD/M)^2)))
    
    ## Assign the distribution type to parameters: 1: normal; 2: lognormal
    rownames(Dat)<-names(pars)
    Data1 <- Dat%>% filter(Dist == 1)
    Data2 <- Dat%>% filter(Dist == 2)
    
    mc1<-list()
    mc2<-list()
    
    for(i in 1:dim(Data1)[1]) {
        mc1[[i]]<-rnormTrunc (N, 
                              min  = qnorm(0.025, mean = Data1$Mean[i], sd = Data1$SD[i]), 
                              max  = qnorm(0.975, mean = Data1$Mean[i], sd = Data1$SD[i]),  
                              mean = Data1$Mean[i], 
                              sd   = Data1$SD[i])
    }
    
    for(i in 1:dim(Data2)[1]) {
        mc2[[i]]<-rlnormTrunc (N, 
                               min = qlnorm(0.025, mean = Data2$logM[i], sd = Data2$logSD[i]), 
                               max = qlnorm(0.975, mean = Data2$logM[i], sd = Data2$logSD[i]), 
                               meanlog = Data2$logM[i], 
                               sdlog   = Data2$logSD[i])
    }
    
    mc1 <-do.call(cbind,mc1)
    mc2 <-do.call(cbind,mc2)
    
    colnames(mc1) <- rownames(Data1)
    colnames(mc2) <- rownames(Data2)
    idata <- cbind(ID = 1:N, mc1, mc2) 
    
    return(as.data.frame(idata))
}



## Source code for Circular bar-plot 

## Loading the package 
#library(tidyverse)

## Function for Circular bar-plot 
Cirplot <- function (data, n = 3) {
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar <- 3
    to_add <- data.frame( matrix(NA, empty_bar*nlevels(factor(data$Chem)), ncol(data)))
    colnames(to_add) <- colnames(data)
    to_add$Chem <- rep(levels(factor(data$Chem)), each=empty_bar)
    data <- rbind(data, to_add)
    data <- data %>% arrange(Chem)
    data$id <- seq(1, nrow(data))
    
    # Get the name and the y position of each label
    label_data <- data
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse(angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data <- data %>% 
        group_by(Chem) %>% 
        summarise(start=min(id), end=max(id) - empty_bar) %>% 
        rowwise() %>% 
        mutate(title=mean(c(start, end)))
    
    # prepare a data frame for grid (scales)
    grid_data <- base_data
    grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1,]
    
    ## Add the font to font database
    windowsFonts("Times" = windowsFont("Times New Roman"))
    
    # Make the plot
    p <- ggplot(data, aes(x=as.factor(id), y=value, fill=Chem)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
        
        geom_bar(aes(x=as.factor(id), y=value, fill=Chem), stat="identity", alpha=0.5) +
        
        # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
        geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        
        # Add text showing the value of each 100/75/50/25 lines
        annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
        
        geom_bar(aes(x=as.factor(id), y=value, fill=Chem), stat="identity", alpha=0.5) +
        ylim(-100,200) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            plot.margin = unit(rep(-1,4), "cm") 
        ) +
        coord_polar() + 
        geom_text(data=label_data, 
                  aes(x=id, y=value+10, label=Pars, hjust=hjust), 
                  color ="black", family = "Times",
                  fontface ="bold",
                  alpha = 1, 
                  size = 6, angle= label_data$angle, inherit.aes = FALSE ) +
        
        # Add base line information
        geom_segment(data=base_data, aes(x = start, y = -2, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
        ## Add drug labels in the circle
        geom_text(data=base_data, aes(x = title, y = -18, label=""), hjust = c(1,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    
    return (p)
}


## Set-up the parameters for plot
windowsFonts(Times=windowsFont("Times New Roman"))

## Plot theme
ptheme<-theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border = element_rect(colour = "black", fill=NA, size=2.5),
    #axis.line.x = element_line(size = 2, linetype = "solid", colour = "black"),
    #axis.line.y = element_line(size = 2, linetype = "solid", colour = "black"),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    ggh4x.axis.ticks.length.minor = rel(0.5),
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),
    axis.ticks.length.x = unit(.25, "cm"),
    axis.ticks = element_line(size = 1),# label of axes
    legend.position='none') 


## Define the plot function
###---------------------------------------------------
MCplot <- function (data, target, Species = "BC", TOL, tdoses, tinterval=24, color="gray", X.limit = NULL, Y.limit=1e-3) {
    
    if (Species == "BC") {
     PDat <- data %>% group_by(Time)%>%
        summarise( Time = Time,
                   Median = switch(target,
                                   "Plasma"  = median(CP),
                                   "Liver"   = median(CL),
                                   "Kidney"  = median(CK),
                                   "Muscle"  = median(CM)
                   ),

                   Ub     = switch(target,
                                   "Plasma"  = quantile(CP, probs  = 0.99),
                                   "Liver"   = quantile(CL, probs  = 0.99),
                                   "Kidney"  = quantile(CK, probs  = 0.99),
                                   "Muscle"  = quantile(CM, probs  = 0.99)
                   ),

                   Lb     = switch(target,
                                   "Plasma"  = quantile(CP, probs  = 0.01),
                                   "Liver"   = quantile(CL, probs  = 0.01),
                                   "Kidney"  = quantile(CK, probs  = 0.01),
                                   "Muscle"  = quantile(CM, probs  = 0.01),
                                   .groups = 'drop')
        )
        
    
    } else {
        
        PDat <- data %>% group_by(Time)%>%
            summarise(Time = Time,
                       Median = switch(target,
                                       "Plasma"  = median(CP),
                                       "Liver"   = median(CL),
                                       "Kidney"  = median(CK),
                                       "Muscle"  = median(CM),
                                       "Milk"    = median(CMilk)
                       ),
                       
                       Ub     = switch(target,
                                       "Plasma"  = quantile(CP, probs  = 0.99),
                                       "Liver"   = quantile(CL, probs  = 0.99),
                                       "Kidney"  = quantile(CK, probs  = 0.99),
                                       "Muscle"  = quantile(CM, probs  = 0.99),
                                       "Milk"    = quantile(CMilk, probs  = 0.99)
                                       
                       ),
                       
                       Lb     = switch(target,
                                       "Plasma"  = quantile(CP, probs  = 0.01),
                                       "Liver"   = quantile(CL, probs  = 0.01),
                                       "Kidney"  = quantile(CK, probs  = 0.01),
                                       "Muscle"  = quantile(CM, probs  = 0.01),
                                       "Milk"    = quantile(CMilk, probs  = 0.01),
                                       .groups = 'drop')
            )
                       
    }
    
    

    A <- PDat %>%filter (Time >= ((tdoses)*tinterval) & round(Ub,3) < TOL)%>% 
                   select(Time)%>%min()
    
    x.intercept <- ifelse (A == "Inf", 0, A)
    
    if (is.null(X.limit) == TRUE) {
    endtime     <- (x.intercept + 30*12*tinterval)/24
    } else {endtime = X.limit}
    
    p <-ggplot (PDat, aes(Time/24, Median))+ 
        geom_ribbon(aes(ymin = Lb, ymax = Ub), fill =color)+
        geom_line  (aes(y = Ub), colour = "black",  lwd = 1,linetype=2) +
        geom_line  (aes(y = Lb), colour = "black",  lwd = 1,linetype=2) +
        geom_line  (lwd=1, colour = "black") +
        scale_x_continuous ( breaks = pretty(0:endtime, n = 5),
                             #minor_breaks = pretty(0:endtime, n=endtime),
                             label  = pretty((0:endtime), n = 5)-((tdoses-1)*tinterval/24),
                             limits = c(0, endtime)) +
        scale_y_log10( 
            labels = function(x) format(x, scientific = TRUE),
            limits = c(Y.limit, NA))+
        annotation_logticks( size = 0.6,
                            sides = "l") + ptheme + labs(x="",y="",face="bold")
    
    p1 <- p + geom_vline(aes(xintercept  = x.intercept/24),
                         size = 1, color = "red", linetype = 2,
                         show.legend = F) 
    
        
    a <- PDat %>% filter (Time > ((tdoses)*tinterval)) %>% select(Ub)
    
    if (max(a$Ub) < TOL | A == "Inf") { p2 = p + geom_line(aes(y = TOL),color = 'black',size = 1, linetype = 'twodash', show.legend = F)}
    else { p2 = p1 + geom_line(aes(y = TOL),color = 'black',size = 1, linetype = 'twodash', show.legend = F)}

   
    return(p2)
}













