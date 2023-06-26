
##---------------------------------------------------------------------------------------------------
# This mrgsolve PBPK model code is developed for the generic model 
# to simulate the concentrations of PFOA, PFOS and PFXS 
# for Beef cattle
#----------------------------------------------------------------------------------------------------

PFOS_Beef_PBPK.code <- '
$PROB @annotated
# Generic PBPK model for PFOA, PFOS and PFXS  
- Author    : Wei-Chun Chou
- Adivisor  : Zhoumeng Lin
- Date      : Initiated in Oct, 2021, finalized in June 2023
- Strucutre : GI tract, Muscle, rest of body, kidneys, liver, venous and arterial blood
- Default physiological parameters value obtained from Lin et al. (2020)

$PARAM @annotated
BW                  : 512     : kg,                  Body weight based on measurement data (Drew et al., 2019)
Htc                 : 0.378   : Unitless             Hematocrit for adult cattle (Lin et al. 2020)
QCC                 : 5.45    : L/h/kg,              Cardiac outputin unanesthetized adult cattle (Lin et al. 2020)
QLCa                : 0.44    : Unitless,            Fractional blood flow to liver (Brown 1997)
QKCa                : 0.11    : Unitless,            Fractional blood flow to kidney (Brown 1997)
QMCa                : 0.28    : Unitless,            Fractional blood flow to muscle (Brown 1997)
QRestCa             : 0.17    : Fraction of blood flow to the rest of body, Value from (total sum equals to 1) (1-QLC-QKC-QFC-QMC)
VLCa                : 0.012   : Unitless,            Fractional liver tissue (Brown 1997)
VKCa                : 0.002   : Unitless,            Fractional kidney tissue (Brown 1997)
VMCa                : 0.36    : L/kg BW,
VRestCa             : 0.587   : Fractional rest of body, Value calculated from 1-VLC-VKC-VFC-VMC-VbloodC
VFilC               : 0.0002  : L/kg BW,             Fraction vol. of filtrate; 10% of Kidney volume; (Worley and Fisher et al., 2015)
VPTCC               : 0.135   : L/kg kidney,         Volume of proximal tubule cells (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)(Hsu et al., 2014)
VbloodCa            : 0.039   : L/kg BW,             Fractional plasma (Davies 1993)
BVKC                : 0.066   : Unitless,            Blood volume fraction of kidney (Brown, 1997)
PL                  : 1.03    : Unitless,            Liver/ plasma PC; (Loccisano et al., 2012) 
PK                  : 0.48    : Unitless,            Kidney/ plasma PC; (Loccisano et al., 2012)
PM                  : 0.08    : Unitless,
PRest               : 0.1526  : Unitless,
Free                : 0.055   : Unitless,
KabsC               : 2.12    : 1/h/kg BW^(-0.25),   Rate of absorption of PFOA in small intestines (Worley et al., 2015)
KurineC             : 0.001   : 1/h/kg BW^(-0.25),
KbileC              : 0.001   : L/h/kg, biliary secretion rate constant for drug,
KehcC               : 0.01    : 1/h/kg, rate constant for the enterohepatic circulation
GFRC                : 0.138   : L/h/kg, Orinral value (177.5 mL/min/m2 body surface area) collected from Murayama et al. 2013; The conversion with (0.1775*60*(0.09*512^0.69))/512
Kdif                : 0.001   : L/h,
KeffluxC            : 0.1     : 1/h/kg BW^(-0.25),
KfecesC             : 0.0137  : 1/h/kg BW^(-0.25),
GE                  : 0.182    : 1/h, Gastric emptying time  (Croubels et al., 2006; Yang et al., 2019)
MW                  : 500.126 : g/mol,               PFOS molecular mass 
Nprotein            : 2e-6    : mg protein/PTC,      Amount of protein in proximal tubule cells (PTCs) (Addis et al., 1936)
RAF_apical          : 0.001   : Unitless,            Relative acitivty factor of apical transpoters, assumed the same as for humans (Chou and Lin, 2019)
Vmax_apical_invitro : 51803   : pmol/mg protein/min, Vmax of apical transporter OATp1a1, initiated based on the value for human (Chou and Lin, 2019) 
Km_apical           : 64.400  : mg/L,                Km of apical transpoter OATp1a1, initiated based on the value for human (Chou and Lin, 2019) 
RAF_baso            : 0.005   : Unitless             Relative activity factor of basolateral transpoters, assumed the same as for humans (Chou and Lin, 2019) 
Vmax_baso_invitro   : 60.99   : pmol/mg protein/min, Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
Km_baso             : 34.500  : mg/L,                Km of basolateral transpoter, averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)


$MAIN
double sumQ        = QLCa + QKCa + QMCa + QRestCa;                    // Sum up cardiac output fraction
double sumV        = VLCa + VKCa + VMCa + VbloodCa + VRestCa;         // Sum up the tissue volumes
double QLC         = QLCa/sumQ;                                       // Adjusted blood flow rate fraction to liver
double QKC         = QKCa/sumQ;                                       // Adjusted blood flow rate fraction to kidney
double QMC         = QMCa/sumQ;                                       // Adjusted blood flow rate fraction to muscle
double QRestC      = 1-QLC-QKC-QMC;                                   // Adjusted blood flow rate fraction to rest of body
double VLC         = VLCa/sumV;                                       // Adjusted fraction of tissue volume of liver
double VKC         = VKCa/sumV;                                       // Adjusted fraction of tissue volume of kidney
double VMC         = VMCa/sumV;                                       // Adjusted fraction of tissue volume of muscle
double VbloodC     = VbloodCa/sumV;                                   // Adjusted fraction of tissue volume of blood
double VRestC      = 1-VLC-VbloodC-VKC-VMC;                           // Adjusted fraction of tissue volume of rest of body
double QC          = QCC*BW;                                          // Cardiac output
double QL          = QLC*QC;                                          // Blood flows in Liver
double QK          = QKC*QC;                                          // Blood flows in Kidney
double QM          = QMC*QC;                                          // Blood flows in Muscle
double QRest       = QRestC*QC;                                       // Blood flows in Rest of body (parent compound submodel)
double VL          = VLC*BW;                                          // Tissues vloume of Liver
double VK          = VKC*BW;                                          // Tissues vloume of Kidney
double VKb         = VK*BVKC;                                
double VM          = VMC*BW;                                          // Tissues vloume of Muscle
double VRest       = VRestC*BW;                                       // Tissues vloume of Rest of body (parent compound submodel)
double Vblood      = VbloodC*BW;                                      // BTissues vloume of lood
double VPlas       = Vblood*(1-Htc);                                  // Tissues vloume of Plasma
double PTC         = VKC*6e7*1000;                           
double VFil        = VFilC*BW;                                
double VPTC        = VK*VPTCC;                               
double Vmax_basoC  = (Vmax_baso_invitro*RAF_baso*PTC*Nprotein*60*(MW/1e12)*1e3);         
double Vmax_apicalC = (Vmax_apical_invitro*RAF_apical*PTC*Nprotein*60*(MW/1e12)*1e3);      
double Vmax_baso   = Vmax_basoC*pow(BW,0.75);         
double Vmax_apical = Vmax_apicalC*pow(BW,0.75);     
double Kurine      = KurineC*pow(BW,(-0.25));             
double Kefflux     = KeffluxC*pow(BW,(-0.25));           
double Kfeces      = KfecesC*pow(BW,(-0.25));
double GFR         = GFRC*512;  //  The body sufrace area equation (0.09 * BW^0.67) in cattle from Berman et al. (2003)                            
double Kabs        = KabsC*pow(BW,(-0.25));
double Kbile       = KbileC*BW;                                       // Biliary elimination rate constant for parent compound
double Kehc        = KehcC*BW;                                        // Rate constant for enterohepatic circulation


// #+ Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE          = (TIME==0)?0:1;  

$CMT A_baso A_apical Adif Aefflux ACL APL APTC AFil Aurine AKb ARest AST ASI
AabsST Afeces AL AM ADOSE AUCCP AUCCL AUCCK AUCCM

$ODE

// Concentrations in the tissues and in the venous palsma leaving each of the tissues
double CPL       = APL/VPlas;                  
double CPL_all   = CPL/Free;                           
double CL        = AL/VL;                                   
double CVL       = CL/PL;                                  
double CKb       = AKb/VKb;                                 
double CVK       = CKb;
double CK        = CVK*PK;
double CM        = AM/VM;
double CVM       = CM/PM;
double CRest     = ARest/VRest;                          
double CVRest    = CRest/PRest;
double Curine    = Aurine;

// Kidney compartment plus 2 subcompartment (Proximal Tubule cells: PTCs, Filtrate: Fil)
// Concentration in kidney, PTCs and fil
double CPTC      = APTC/VPTC;                             
double CFil      = AFil/VFil;                             

// Virtural compartment: 
// Basolateral (baso) 
// transport, Diffusion (dif), 
// Apical (apical) transport, and 
// efflux clearance (efflux)
// Clerance (CL) via glormerular filtration 

double RA_baso        = (Vmax_baso*CKb)/(Km_baso + CKb);                                             
double RA_apical      = (Vmax_apical*CFil)/(Km_apical + CFil);                                    
double Rdif           = Kdif*(CKb - CPTC);                                                            
double RA_efflux      = Kefflux*APTC;                                                             
double RCL            = CPL_all*GFR*Free;                                                                   
double RPTC           = Rdif + RA_apical + RA_baso - RA_efflux;                                         
double Rurine         = Kurine*AFil;                                                                
double RFil           = RCL - RA_apical - Rurine;                                        
double RKb            = QK*(CPL_all - CVK)*Free - CPL_all*GFR*Free - Rdif - RA_baso;                               
double Rrest          = QRest*(CPL_all - CVRest)*Free;                                                      
double RST            = -GE*AST;                                                             
double RSI            =  GE*AST - Kabs*ASI - Kfeces*ASI- Kehc*ASI + Kbile*AL;                                                
double RabsST         = Kabs*ASI;                                                                     
double Rfeces         = Kfeces*ASI;                                                      
double Rbile          = Kbile*AL;   
double RM             = QM*(CPL_all - CVM)*Free;
double RL             = QL*(CPL_all - CVL)*Free + RabsST - Rbile + Kehc*ASI;                                 
double RPL            = (QL*CVL + QK*CVK + QM*CVM + QRest*CVRest)*Free - QC*CPL_all*Free + RA_efflux;  


// ## ODE
dxdt_A_baso           = RA_baso;                                                                      
dxdt_A_apical         = RA_apical;                                                                  
dxdt_Adif             = Rdif;                                                                           
dxdt_Aefflux          = RA_efflux;                                                                    
dxdt_ACL              = RCL;                                                                             
dxdt_APL              = RPL;                                                              
dxdt_APTC             = RPTC;                                                                           
dxdt_AFil             = RFil;                                                                           
dxdt_Aurine           = Rurine;                                                                       
dxdt_AKb              = RKb;                                                                             
dxdt_ARest            = Rrest;                                                                         
dxdt_AST              = RST;
dxdt_ASI              = RSI; 
dxdt_AabsST           = RabsST;                                                                       
dxdt_Afeces           = Rfeces;                                                                       
dxdt_AL               = RL;                                                                               
dxdt_AM               = RM;
dxdt_AUCCP            = CPL_all;
dxdt_AUCCL            = CL;
dxdt_AUCCK            = CK;
dxdt_AUCCM            = CM;



// #+ Virtural compartment; input dose for mass balance
dxdt_ADOSE            = 0;

// Mass balance equations
double Qbal           = QC - QL - QK - QM - QRest;
double Vbal           = BW - Vblood - VL - VM - VRest - VK;
double Tmass          = APL + AL + AM + ARest + AKb + AFil + APTC + AST + ASI;
double Loss           = Aurine + Afeces;
double ATotal         = Tmass + Loss;
double Mbal           = ADOSE - ATotal*KDOSE; 

$TABLE
capture Plasma        = CPL_all;
capture Liver         = CL;
capture Kidney        = CK;
capture Muscle        = CM;
capture Urine         = Aurine;
capture AUC_CP        = AUCCP;
capture AUC_CL        = AUCCL;
capture AUC_CK        = AUCCK;
capture AUC_CM        = AUCCM;
//capture MBQ          = Qbal;
//capture MBV          = Vbal;
'
