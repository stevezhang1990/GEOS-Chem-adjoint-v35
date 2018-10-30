# PHOTOLYSIS REACTIONS - MASTER RATEFILE - Paul Brown, Oliver Wild & David Rowley
# Centre for Atmospheric Science, Cambridge, U.K.  Release date:  22 November 1993
# SCCS version information: @(#)photol.d	1.2 5/11/94
# 
# Modified for Harvard chemistry: several reactions added, re-ordered per chem.dat
# Also putting in the Harvard names in col 1, the UCI x-sec names in last col !!!
# 					-Prashant Murti [4/13/98]
#
# The new peroxide recycling now activates the following photolysis species:
#   GP,IAP,INPN,ISN1,ISNP,MAOP,MRP,PP,PRPN,RIP,VRP.
# Also be sure to set parameter JPMAX = 55 in "cmn_fj.h". 
#                                       - Randall Martin & Bob Yantosca [12/20/00]
# New updates from FASTJX.(jmao,ccarouge, 04/20/09)
#
# Harvard species             Products - UCI notation                              UCI xsec
# ===============           ===============================                        ========
    1 H2O        PHOTON     OH         HO2                  0.00E+00  0.00      0.0         
    2 HO2        PHOTON     OH         O(3P)                0.00E+00  0.00      0.0         
    3 O2         PHOTON     O(3P)      O(3P)                0.00E+00  0.00    100.0  O2     
    4 O3_P       PHOTON     O2         O(3P)                0.00E+00  0.00    100.0  O3     
    5 O3         PHOTON     O2         O(1D)                0.00E+00  0.00    100.0  O3_1d  
    6 NO2        PHOTON     NO         O(3P)                0.00E+00  0.00    100.0  NO2    
    7 H2O2       PHOTON     OH         OH                   0.00E+00  0.00    100.0  H2O2   
    8 MP         PHOTON     HCHO       OH         HO2       0.00E+00  0.00    100.0  ROOH   
    9 CH2O       PHOTON     CO         HO2        HO2       0.00E+00  0.00    100.0  HCHO=H+
   10 CH2O       PHOTON     CO         H2                   0.00E+00  0.00    100.0  HCHO=H2
   11 HNO3       PHOTON     OH         NO2                  0.00E+00  0.00    100.0  HONO2  
   12 HNO2       PHOTON     OH         NO                   0.00E+00  0.00    100.0  HONO   
   13 HNO4       PHOTON     OH         NO3                  0.00E+00  0.00     33.3  HO2NO2 
   14 HNO4       PHOTON     HO2        NO2                  0.00E+00  0.00     66.7  HO2NO2 
   15 NO3        PHOTON     NO         O2                   0.00E+00  0.00    100.0  NO3=O2+
   16 NO3        PHOTON     NO2        O(3P)                0.00E+00  0.00    100.0  NO3=O+ 
   17 N2O5       PHOTON     NO3        NO         O(3P)     0.00E+00  0.00      0.0  N2O5   
   18 N2O5       PHOTON     NO3        NO2                  0.00E+00  0.00    100.0  N2O5   
   19 ALD2       PHOTON     CH4        CO                   0.00E+00  0.00    100.0  Acet=R+
   20 ALD2       PHOTON     MeOO       HO2        CO        0.00E+00  0.00    100.0  Acet=RO
   21 PAN        PHOTON     MeCO3      NO2                  0.00E+00  0.00    100.0  PAN    
   22 RCHO       PHOTON     EtO2       HO2        CO        0.00E+00  0.00    100.0  RCHO   
   23 ACET       PHOTON     MeCO3      MeOO                 0.00E+00  0.00    100.0  AcetA 
   24 ACET       PHOTON     MeOO       MeOO       CO        0.00E+00  0.00    100.0  AcetB 
   25 MEK        PHOTON     MeCO3      EtOO                 0.00E+00  0.00    100.0  EtCOMe
   26 MNO3       PHOTON     HCHO       H2O        NO2       0.00E+00  0.00    100.0  MeNO3  
   27 GLYC       PHOTON     HCHO       HO2        CO        0.00E+00  0.00    100.0  HOMeCHO
   28 GLYX       PHOTON     H2         CO         HCHO      0.00E+00  0.00    100.0  Glyxla 
   29 GLYX       PHOTON     CO         HO2                  0.00E+00  0.00    100.0  Glyxlb 
   30 MGLY       PHOTON     MeCO3      CO         HO2       0.00E+00  0.00    100.0  MeCOCHO
   31 MGLY       PHOTON     Acet       CO                   0.00E+00  0.00      0.0  MeCOCHO
   32 MVK        PHOTON     PRPE       CO                   0.00E+00  0.00     60.0  MeCOVi 
   33 MVK        PHOTON     MeCO3      HCHO       CO  HO2   0.00E+00  0.00     20.0  MeCOVi 
   34 MVK        PHOTON     MeOO       MAO3                 0.00E+00  0.00     20.0  MeCOVi 
   35 MACR       PHOTON     MAO3       HO2                  0.00E+00  0.00     50.0  MACR   
   36 MACR       PHOTON     CO  HO2  MGLY  HO2  MeCO3  HCHO 0.00E+00  0.00     50.0  MACR   
   37 HAC        PHOTON     MeCO3      HCHO       HO2       0.00E+00  0.00    100.0  AcetA  
   38 ETP        PHOTON     OH         HO2        Acet      0.00E+00  0.00    100.0  ROOH   
   39 RA3P       PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH   
   40 RB3P       PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH   
   41 R4P        PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH   
   42 RP         PHOTON     OH         HO2        Acet      0.00E+00  0.00    100.0  ROOH   
   43 R4N2       PHOTON     NO2 MeCOMe MEK MO2 HO2 ALD2 ... 0.00E+00  0.00    100.0  MeNO3  
   44 MAP        PHOTON     OH         MO2                  0.00E+00  0.00    100.0  ROOH   
   45 INPN       PHOTON     OH         HO2        RCHO NO2  0.00E+00  0.00    100.0  ROOH
   46 PRPN       PHOTON     OH         HO2        RCHO NO2  0.00E+00  0.00    100.0  ROOH
   47 PP         PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH
   48 GP         PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH
   49 GLP        PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH
   50 RIP        PHOTON     OH HO2 CH2O MVK MACR RIO1 IALD  0.00E+00  0.00    100.0  ROOH
   51 IAP        PHOTON     OH HO2  CO  H2  HAC  GLYC  MGLY 0.00E+00  0.00    100.0  ROOH
   52 ISNP       PHOTON     OH         HO2        RCHO NO2  0.00E+00  0.00    100.0  ROOH
   53 VRP        PHOTON     OH  HO2  CH2O  MCO3  GLYC  MGLY 0.00E+00  0.00    100.0  ROOH
   54 MRP        PHOTON     OH  HO2  MGLY  HAC  CO  CH2O    0.00E+00  0.00    100.0  ROOH 
   55 MAOP       PHOTON     OH         HO2        RCHO      0.00E+00  0.00    100.0  ROOH
 9999                                                       0.00E-00  0.00      0.0         



 NOTES:                                                                          
 -----                                                                           
[4/15/98]
Oliver Wild: All reaction data from JPL '97, IUPAC IV. IUPAC V is soon
expected. - ppm


  All reaction data taken from IUPAC supplement IV unless otherwise indicated.

  JPL - data from JPL (latest assessment as far as possible)
                                                                                 
  ? - reaction products unknown
  * - user strongly advised to consult source material
  B - branching ratio assumed equal for all channels in the absence of more information
  U - upper limit for rate coefficient


 Changes since 08/3/93 release:
  O now written as O(3P)
                                                                                 
(Note that the second of the acetaldehyde channels above occurs at wavelengths
less than 289 nm, and therefore doesn't appear in the Fast-J region at all - 
I've simply included it here for completeness) - [from Oliver, 3/7/98]
