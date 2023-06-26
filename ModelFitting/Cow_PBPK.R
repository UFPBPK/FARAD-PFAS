##---------------------------------------------------------------------------------------------------
# This mrgsolve PBPK model code is developed for the generic model 
# to simulate the concentrations of PFOA, PFOS and PFXS 
# for dairy cow
#----------------------------------------------------------------------------------------------------

PFOS_Cow_PBPK.code <- '
$PROB @annotated
# Generic PBPK model for PFOA, PFOS and PFXS  
- Author    : Wei-Chun Chou
- Adivisor  : Zhoumeng Lin
- Date      : June, 2023
- Strucutre : GI tract, Muscle, rest of body, fat, kidneys, liver, venous and arterial blood
- Default physiological parameters value obtained from Li et al. 2017 and Lin et al. (2020)

$PARAM @annotated
BW                  : 583     : kg,                  Body weight based on measurement data (Kowalczyk et al., 2013)
Htc                 : 0.29    : Unitless,            Hematocrit for Dairy cow (Herman et al. 2018)
QCC                 : 5.45    : L/h/kg,              Cardiac outputin unanesthetized adult cattle (Lin et al. 2020)
QLCa                : 0.58    : Unitless,            Fraction of blood flow of liver tissue
QKCa                : 0.07    : Unitless,            Fraction of blood flow of kidney tissue 
QMCa                : 0.18    : Unitless,            Fraciton of blood flow of muscle tissue (Li et al. 2017) 
QUCa                : 0.13    : Unitless,            Fraciton of blood flow of udder tissue 
QRestCa             : 0.04    : Fraction of blood flow to the rest of body, Value from (total sum equals to 1) (1-QLCa-QKCa-QMCa-QUCa)
VLCa                : 0.0135  : Unitless,            Fraction of tissue voloume of liver (dairy cows in breeds of Jersey)
VKCa                : 0.0022  : Unitless,            Fraction of tissue voloume of kidney (dairy cows in breeds of Jersey)
VMCa                : 0.27    : L/kg BW,             Fraction of tissue voloume of muscle (Li et al. 2017) 
VUCa                : 0.0314  : L/kg BW,             Fraction of tissue voloume of udder (dairy cows in breeds of Jersey)
VbloodCa            : 0.0433  : L/kg BW,             Fraction of tissue voloume of blood (dairy cows in breeds of Jersey)
VRestCa             : 0.6396  : Fractional rest of body, Value calculated from 1-VLCa-VKCa-VMCa-VUCa-VbloodCa
BVKC                : 0.066   : Unitless,            Blood volume fraction of kidney (Brown, 1997)
VFilC               : 0.0002  : L/kg BW,             Fraction vol. of filtrate; 10% of Kidney volume; (Worley and Fisher et al., 2015)
VPTCC               : 0.135   : L/kg kidney,         Volume of proximal tubule cells (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)(Hsu et al., 2014)
PL                  : 1.6534  : Unitless,            
PK                  : 1.35    : Unitless,            
PM                  : 0.0772  : Unitless,
PU                  : 0.2     : Unitless,
PRest               : 0.05    : Unitless,
Free                : 0.0152  : Unitless,
KabsC               : 2.12    : 1/h/kg BW^(-0.25),
KurineC             : 0.001   : 1/h/kg BW^(-0.25),
KbileC              : 0.001   : L/h/kg, biliary secretion rate constant for drug,
KehcC               : 0.01    : 1/h/kg, rate constant for the enterohepatic circulation
GFRC                : 0.138   : L/h/kg, Orinral value (177.5 mL/min/m2 body surface area) collected from Murayama et al. 2013; The conversion with (0.1775*60*(0.09*512^0.69))/512
Kdif                : 0.001   : L/h,
KeffluxC            : 5       : 1/h/kg BW^(-0.25),
KfecesC             : 0.002   : 1/h/kg BW^(-0.25),
GE                  : 0.182   : 1/h, Gastric emptying time  (Croubels et al., 2006; Yang et al., 2019)
MW                  : 500.126 : g/mol,               Molecular mass-PFOS MW as default value 
Nprotein            : 2e-6    : mg protein/PTC,      Amount of protein in proximal tubule cells (PTCs) (Addis et al., 1936)
RAF_apical          : 0.001   : Unitless,            Relative acitivty factor of apical transpoters, assumed the same as for humans (Chou and Lin, 2019)
Vmax_apical_invitro : 51803   : pmol/mg protein/min, Vmax of apical transporter OATp1a1, initiated based on the value for human (Chou and Lin, 2019) 
Km_apical           : 64.400  : mg/L,                Km of apical transpoter OATp1a1, initiated based on the value for human (Chou and Lin, 2019) 
RAF_baso            : 0.005   : Unitless             Relative activity factor of basolateral transpoters, assumed the same as for humans (Chou and Lin, 2019) 
Vmax_baso_invitro   : 60.99   : pmol/mg protein/min, Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
Km_baso             : 34.500  : mg/L,                Km of basolateral transpoter, averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOS (Worley and Fisher, 2015)
VMilksp             : 18.72   : L,
PMilkM              : 0.9587  : Unitless,
QPAMilkC            : 0.5     : L/h/kg BW^0.75,
Kmilking            : 0.1801  : L/h,

$MAIN
double sumQ        = QLCa + QKCa + QMCa + QUCa + QRestCa;             // Sum up cardiac output fraction
double sumV        = VLCa + VKCa + VMCa + VUCa + VbloodCa + VRestCa;         // Sum up the tissue volumes
double QLC         = QLCa/sumQ;                                       // Adjusted blood flow rate fraction to liver
double QKC         = QKCa/sumQ;                                       // Adjusted blood flow rate fraction to kidney
double QMC         = QMCa/sumQ;                                       // Adjusted blood flow rate fraction to muscle
double QUC         = QUCa/sumQ;                                       // Adjusted blood flow rate fraction to muscle
double QRestC      = 1-QLC-QKC-QUC-QMC;                                // Adjusted blood flow rate fraction to rest of body
double VLC         = VLCa/sumV;                                       // Adjusted fraction of tissue volume of liver
double VKC         = VKCa/sumV;                                       // Adjusted fraction of tissue volume of kidney
double VMC         = VMCa/sumV;                                       // Adjusted fraction of tissue volume of muscle
double VUC         = VUCa/sumV;                                       // Adjusted fraction of tissue volume of muscle
double VbloodC     = VbloodCa/sumV;                                   // Adjusted fraction of tissue volume of blood
double VRestC      = 1-VLC-VbloodC-VKC-VMC-VUC;                       // Adjusted fraction of tissue volume of rest of body

double QC          = QCC*BW;               
double QK          = QKC*QC;                                  
double QL          = QLC*QC;                                 
double QM          = QMC*QC;
double QU          = QUC*QC;
double QPAMilk     = QPAMilkC*pow(BW,0.75);
double QRest       = QC*QRestC;     

double VPL         = VbloodC*BW;
double VL          = VLC*BW;                                  
double VK          = VKC*BW;                                   
double VM          = VMC*BW;
double VU          = VUC*BW;
double VFil        = VFilC*BW;                                
double VPTC        = VK*VPTCC;                               
double VKb         = VK*BVKC;                                
double VRest       = VRestC*BW;;      

double PTC         = VKC*6e7*1000;                           
double Vmax_basoC  = (Vmax_baso_invitro*RAF_baso*PTC*Nprotein*60*(MW/1e12)*1e3);         
double Vmax_apicalC = (Vmax_apical_invitro*RAF_apical*PTC*Nprotein*60*(MW/1e12)*1e3);      
double Vmax_baso   = Vmax_basoC*pow(BW,0.75);         
double Vmax_apical = Vmax_apicalC*pow(BW,0.75);     
double Kurine      = KurineC*pow(BW,(-0.25));             
double Kefflux     = KeffluxC*pow(BW,(-0.25));           
double Kfeces      = KfecesC*pow(BW,(-0.25));
double GFR         = GFRC*BW;  //  The body sufrace area equation (0.09 * BW^0.67) in cattle from Berman et al. (2003)                            
double Kabs        = KabsC*pow(BW,(-0.25));                 
double Kbile       = KbileC*BW;                                       // Biliary elimination rate constant for parent compound
double Kehc        = KehcC*BW;                                        // Rate constant for enterohepatic circulation

// #+ Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE          = (TIME==0)?0:1;  

$CMT A_baso A_apical Adif Aefflux ACL APL APTC  AFil  Aurine AKb  ARest  AST ASI
AabsST Afeces AL AMilk AMilking AM AU ADOSE AUCCP AUCCL AUCCK AUCCM AUCCMilk

$ODE

// Concentrations in the tissues and in the venous palsma leaving each of the tissues
double CPL       = APL/VPL;                  
double CPL_all   = CPL/Free;                           
double CL        = AL/VL;                                   
double CVL       = CL/PL;                                  
double CKb       = AKb/VKb;                                 
double CVK       = CKb;
double CK        = CVK*PK;
double CM        = AM/VM;
double CVM       = CM/PM;
double CU        = AU/VU;
double CVU       = CU/PU;
double CRest     = ARest/VRest;                          
double CVRest    = CRest/PRest;
double CMilk     = AMilk/VMilksp;

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
double RFil           = CPL_all*GFR*Free - RA_apical - Rurine;                                        
double RKb            = QK*(CPL_all - CVK)*Free - CPL_all*GFR*Free - Rdif - RA_baso;                               
double Rrest          = QRest*(CPL_all - CVRest)*Free;                                                      
double RST            = -GE*AST;                                                             
double RSI            =  GE*AST - Kabs*ASI - Kfeces*ASI- Kehc*ASI + Kbile*AL;                                                
double RabsST         = Kabs*ASI;                                                                     
double Rfeces         = Kfeces*ASI;
double Rbile          = Kbile*AL;   
double RM             = QM*(CPL_all - CVM)*Free;
double RU             = QU*(CPL_all - CVU)*Free - QPAMilk*(CVU*Free - CMilk/PMilkM);
double RL             = QL*(CPL_all - CVL)*Free + RabsST - Rbile + Kehc*ASI;                                 
double RPL            = (QL*CVL + QK*CVK + QM*CVM + QU*CVU + QRest*CVRest)*Free - QC*CPL_all*Free + RA_efflux;  
double RMilking       = Kmilking*CMilk;
double RMilk          = QPAMilk*(CVU*Free - CMilk/PMilkM) - RMilking;

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
dxdt_AMilk            = RMilk;
dxdt_AMilking         = RMilking;
dxdt_AM               = RM;
dxdt_AU               = RU;
dxdt_AUCCP            = CPL_all;
dxdt_AUCCL            = CL;
dxdt_AUCCK            = CK;
dxdt_AUCCM            = CM;
dxdt_AUCCMilk         = CMilk;

// #+ Virtural compartment; input dose for mass balance
dxdt_ADOSE            = 0;

// Mass balance equations
double Qbal           = QC - QL - QK - QM - QU - QRest;
double Vbal           = BW - VPL - VL - VM - VU - VRest - VK;
double Tmass          = APL + AL + AM + AU + ARest + AKb + AFil + APTC + AST + AMilk;
double Loss           = Aurine + Afeces + AMilking;
double ATotal         = Tmass + Loss;
double Mbal           = ADOSE - ATotal*KDOSE; 

$TABLE
capture Plasma        = CPL/Free;
capture Milk          = CMilk;
capture Liver         = CL;
capture Kidney        = CK;
capture Muscle        = CM;
capture Urine         = Aurine;
capture AUC_CP        = AUCCP;
capture AUC_CL        = AUCCL;
capture AUC_CK        = AUCCK;
capture AUC_CM        = AUCCM;
capture AUC_CMilk     = AUCCMilk;
'