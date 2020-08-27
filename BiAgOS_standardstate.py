#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:17:41 2018

@author: laurenwalters
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import random
#For saving/importing data
from numpy import asarray
from numpy import save
from numpy import load

#Created by Lauren Walters, 2018-2020
#Contributions by Liang Feng Huang
#For reactions in aqueous conditions


#find out how much detail you want in your graph
#n=input("Enter the mesh grid detail you want, suggested (30-140): ")
n=30;

#Constants
R=8.31447;         #kJ/(mol*K)
T=298.15;          #K
F= 9.648533*10**4; #kJ/(V*mol)
P=1;               #bar, 10^5*Pa
eta=6
nI=10**-eta;         #Activity Concentration

#Array showing the composition of Ag:Bi:S
composition=np.array([1,1,1])

#pH Range and Constants
lowpH = -2;
highpH = 16;
pHrange = int;
pHrange = highpH-lowpH;
pHcount = pHrange/n; #used to iterate through pH range

#Applied Potential Range and Constants
Ulow = -1.5;          #V
Uhigh = 1.5;          #V
Urange  = Uhigh-Ulow; #V
Ucount = Urange/n;    #used to iterate through U (energy) range


###############################################################################
######################## DFT CALCULATIONS #####################################
###############################################################################
#Electronic Energies in eV/f.u.
#PBEsol with SOC
Ee_Bi= -5.114928333;    
Ee_Bi2O3= -31.163316;   
Ee_Bi2O5= -40.1344765;   
Ee_Bi2O4=-36.7221975;      
Ee_Bi4O7=-68.40888;      
#PBEsol
Ee_Ag=-3.3376147;          
Ee_AgO=-8.8090755;     
Ee_Ag2O=-12.03654;      
Ee_Ag2O3=-21.5480535;
#PBEsol
Ee_O2=-10.281123
Ee_S= -4.391811875;
###############################################################################
########### MULTICOMPONENT SPECIES ############################################
#Calculated with PBEsol+SOC
Ee_Ag2S=-11.127419;
Ee_Ag8S=-29.3975915;

Ee_Ag5BiO4=-46.0621275;
Ee_Ag2BiO3=-30.539865;
Ee_AgBiO3=-26.74457833;

Ee_BiS2=-14.6172585;
Ee_Bi2S3=-24.878388;
Ee_Bi2S2O=-27.2327565;
Ee_Bi2SO4_3=-109.35902;
Ee_Bi14OS24=-247.57619;
Ee_Bi2SO2=-29.50652;

Ee_AgBiS2=-18.05386933;
Ee_AgBi3S5=-42.968627;
Ee_AgBi6S9=-77.70144;
Ee_AgBiSO=-20.1930165;
###############################################################################
###### Vibrational Energy #####################################################
###############################################################################
#Vibrational Energies in eV/f.u.
#From PBEsol Phonon Calculations
Fvib_O2=-0.272; 
F_rot_trans_O2=0.099
Ftot_O2=Fvib_O2+F_rot_trans_O2
F_H = .202

Fvib_S=-0.0091266451372

Fvib_AgO=  0.01360462162 
Fvib_Ag2O=  -0.07478725194
Fvib_Ag2O3= 0.17837389128
Fvib_Ag= -0.038443004392

Fvib_Bi=-0.0761976993239
Fvib_Bi2O3=-0.057653546889
Fvib_Bi2O5=0.14677315404
Fvib_Bi2O4=0.12231438709
Fvib_Bi4O7=0.08741679245

Fvib_Ag2S=-0.1186037340461
Fvib_Ag8S=-0.4056447459661

Fvib_Ag5BiO4=-0.12861399053
Fvib_Ag2BiO3=-0.098643737028
Fvib_AgBiO3=0.0140751179482

Fvib_BiS2=-0.063943629448
Fvib_Bi2S3=-0.1428187610337
Fvib_Bi2S2O=-0.08193190191
Fvib_Bi2SO4_3=0.81266278392
Fvib_Bi14OS24=0.02990373431
Fvib_Bi2SO2=-0.0265520338422

Fvib_AgBiS2=-0.15368879695
Fvib_AgBi3S5=-0.40670413783
Fvib_AgBi6S9=-0.58247934744
Fvib_AgBiSO=-0.0447490075112
###############################################################################
###   Compounds-Calculate the formation energies   ############################
###############################################################################
#Free Energies of Formation in eV/f.u.
dGf_AgO= ((Ee_AgO)+Fvib_AgO)       -(Ee_Ag+Fvib_Ag)-0.5*((Ee_O2)-Ftot_O2); 
dGf_Ag2O= ((Ee_Ag2O)+Fvib_Ag2O)    -2.0*(Ee_Ag+Fvib_Ag)-0.5*((Ee_O2)-Ftot_O2); 
dGf_Ag2O3= ((Ee_Ag2O3)+Fvib_Ag2O3) -2.0*(Ee_Ag+Fvib_Ag)-1.5*((Ee_O2)-Ftot_O2);

dGf_Bi2O3= ((Ee_Bi2O3)+Fvib_Bi2O3) -2.0*(Ee_Bi+Fvib_Bi)-1.5*(Ee_O2-Ftot_O2); 
dGf_Bi2O5= ((Ee_Bi2O5)+Fvib_Bi2O5) -2.0*(Ee_Bi+Fvib_Bi)-2.5*(Ee_O2-Ftot_O2);
dGf_Bi2O4= ((Ee_Bi2O4)+Fvib_Bi2O4) -2.0*(Ee_Bi+Fvib_Bi)-2.0*(Ee_O2-Ftot_O2);
dGf_Bi4O7= ((Ee_Bi4O7)+Fvib_Bi4O7) -4.0*(Ee_Bi+Fvib_Bi)-3.5*(Ee_O2-Ftot_O2);

dGf_Ag2S=(Ee_Ag2S+Fvib_Ag2S)       -2*(Ee_Ag+Fvib_Ag)-(Ee_S+Fvib_S);
dGf_Ag8S=(Ee_Ag8S+Fvib_Ag8S)       -8*(Ee_Ag+Fvib_Ag)-(Ee_S+Fvib_S);

dGf_Ag5BiO4=(Ee_Ag5BiO4+Fvib_Ag5BiO4) -5*(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-2.0*((Ee_O2)-Ftot_O2);
dGf_Ag2BiO3=(Ee_Ag2BiO3+Fvib_Ag2BiO3) -2*(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-1.5*((Ee_O2)-Ftot_O2);
dGf_AgBiO3=(Ee_AgBiO3+Fvib_AgBiO3)    -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-1.5*((Ee_O2)-Ftot_O2);

dGf_BiS2=(Ee_BiS2+Fvib_BiS2)                  -(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S);
dGf_Bi2S3=(Ee_Bi2S3+Fvib_Bi2S3)               -2*(Ee_Bi+Fvib_Bi)-3*(Ee_S+Fvib_S);
dGf_Bi2S2O=(Ee_Bi2S2O+Fvib_Bi2S2O)            -2*(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);
dGf_Bi2SO4_3=(Ee_Bi2SO4_3+Fvib_Bi2SO4_3)      -2*(Ee_Bi+Fvib_Bi)-3*(Ee_S+Fvib_S)-6.0*((Ee_O2)-Ftot_O2);
dGf_Bi14OS24=(Ee_Bi14OS24+Fvib_Bi14OS24)      -14*(Ee_Bi+Fvib_Bi)-24*(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);
dGf_Bi2SO2=(Ee_Bi2SO2+Fvib_Bi2SO2)            -2*(Ee_Bi+Fvib_Bi)-(Ee_S+Fvib_S)-1.0*((Ee_O2)-Ftot_O2);

dGf_AgBiS2=(Ee_AgBiS2+Fvib_AgBiS2)     -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S);
dGf_AgBi3S5=(Ee_AgBi3S5+Fvib_AgBi3S5)  -(Ee_Ag+Fvib_Ag)-3*(Ee_Bi+Fvib_Bi)-5*(Ee_S+Fvib_S);
dGf_AgBi6S9=(Ee_AgBi6S9+Fvib_AgBi6S9)  -(Ee_Ag+Fvib_Ag)-6*(Ee_Bi+Fvib_Bi)-9*(Ee_S+Fvib_S);
dGf_AgBiSO=(Ee_AgBiSO+Fvib_AgBiSO)     -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);

#Set the reference values
dGf_Ag=0.0;
dGf_Bi=0.0;
dGf_S=0.0;
###############################################################################
###############################################################################
###############################################################################


###############################################################################
############## Aqueous Ion Free Energies of Formation #########################
#Free Energies of Formation in eV/f.u.
##Elemental Silver Species
dGf_Ag_Plus=     0.799201
dGf_Ag_2Plus=    2.7796
dGf_AgO_Plus=    2.33733
dGf_AgO_Minus=  -0.238069
dGf_AgOH=       -0.9535
dGf_AgOH2_Minus=-2.696784

##Elemental Bismuth Species
dGf_Bi_3Plus= 0.6430898
dGf_BiOH_2Plus=  -1.6968378
dGf_BiO_Plus=  -1.4977965

#Elemental Sulphur Species
dGf_H2S=-0.283601
dGf_HS_Minus=0.13053
dGf_S_2Minus=0.9521892
dGf_S2_2Minus=0.8563979
dGf_S3_2Minus=0.7791664
dGf_S4_2Minus=0.7204948
dGf_S5_2Minus=0.6803396
dGf_H2S2O3=-5.6329986
dGf_HS2O3_Minus=-5.6156529
dGf_S2O3_2Minus=-5.515915
dGf_S5O6_2Minus=-9.9087
dGf_S4O6_2Minus=-10.5939
dGf_HS2O4_Minus=-6.13203282
dGf_S2O4_2Minus=-5.9842
dGf_S3O6_2Minus=-9.930382
dGf_H2SO3=-5.580528
dGf_HSO3_Minus=-5.464
dGf_SO3_2Minus=-5.03457
dGf_S2O6_2Minus=-10.02
dGf_H2SO4=-7.6901922
dGf_HSO4_Minus=-7.8029389
dGf_SO4_2Minus=-7.6901922
dGf_S2O8_2Minus=-11.361
dGf_HSO5_Minus=     -6.60739025
dGf_S2O5_2Minus=    -8.195817793

#Water
dGf_H2O=-2.458;
###############################################################################
###############################################################################
###############################################################################

###############################################################################
#############   CONVERT from eV to kJ/mol  ####################################
###############################################################################
dGf_Ag=         dGf_Ag*F;   
dGf_AgO=        dGf_AgO*F; 
dGf_Ag2O=       dGf_Ag2O*F;        
dGf_Ag2O3=      dGf_Ag2O3*F;
dGf_Ag_Plus=    dGf_Ag_Plus*F;
dGf_Ag_2Plus=   dGf_Ag_2Plus*F;
dGf_AgO_Plus=   dGf_AgO_Plus*F;
dGf_AgO_Minus=  dGf_AgO_Minus*F;
dGf_AgOH=       dGf_AgOH*F;
dGf_AgOH2_Minus=dGf_AgOH2_Minus*F;

dGf_Bi=   dGf_Bi*F; 
dGf_Bi2O3=  dGf_Bi2O3*F; 
dGf_Bi2O5= dGf_Bi2O5*F; 
dGf_Bi2O4=dGf_Bi2O4*F;
dGf_Bi4O7=dGf_Bi4O7*F;
dGf_Bi_3Plus= dGf_Bi_3Plus*F;
dGf_BiOH_2Plus=  dGf_BiOH_2Plus*F;
dGf_BiO_Plus=  dGf_BiO_Plus*F;

dGf_S=   dGf_S*F;               #eV/unit
dGf_H2S=dGf_H2S*F;               #eV/unit
dGf_HS_Minus=dGf_HS_Minus*F;           #eV/unit
dGf_S_2Minus=dGf_S_2Minus*F;         #eV/unit
dGf_S2_2Minus=dGf_S2_2Minus*F;        #eV/unit
dGf_S3_2Minus=dGf_S3_2Minus*F;        #eV/unit
dGf_S4_2Minus=dGf_S4_2Minus*F;        #eV/unit
dGf_S5_2Minus=dGf_S5_2Minus*F;        #eV/unit
dGf_H2S2O3=dGf_H2S2O3*F;          #eV/unit
dGf_HS2O3_Minus=dGf_HS2O3_Minus*F;     #eV/unit
dGf_S2O3_2Minus=dGf_S2O3_2Minus*F;      #eV/unit
dGf_S5O6_2Minus=dGf_S5O6_2Minus*F;        #eV/unit
dGf_S4O6_2Minus=dGf_S4O6_2Minus*F;       #eV/unit
dGf_HS2O4_Minus=dGf_HS2O4_Minus*F;    #eV/unit
dGf_S2O4_2Minus=dGf_S2O4_2Minus*F;        #eV/unit
dGf_S3O6_2Minus=dGf_S3O6_2Minus*F;      #eV/unit
dGf_H2SO3=dGf_H2SO3*F;            #eV/unit
dGf_HSO3_Minus=dGf_HSO3_Minus*F;          #eV/unit
dGf_SO3_2Minus=dGf_SO3_2Minus*F;        #eV/unit
dGf_S2O6_2Minus=dGf_S2O6_2Minus*F;         #eV/unit
dGf_H2SO4=dGf_H2SO4*F;           #eV/unit
dGf_HSO4_Minus=dGf_HSO4_Minus*F;      #eV/unit
dGf_SO4_2Minus=dGf_SO4_2Minus*F;      #eV/unit
dGf_S2O8_2Minus=dGf_S2O8_2Minus*F;        #eV/unit
dGf_HSO5_Minus=dGf_HSO5_Minus*F;
dGf_S2O5_2Minus=dGf_S2O5_2Minus*F;

dGf_Ag2S=dGf_Ag2S*F;
dGf_Ag8S=dGf_Ag8S*F;

dGf_Ag5BiO4=dGf_Ag5BiO4*F;
dGf_Ag2BiO3=dGf_Ag2BiO3*F;
dGf_AgBiO3=dGf_AgBiO3*F;

dGf_BiS2=dGf_BiS2*F;
dGf_Bi2S3=dGf_Bi2S3*F;
dGf_Bi2S2O=dGf_Bi2S2O*F;
dGf_Bi2SO4_3=dGf_Bi2SO4_3*F;
dGf_Bi14OS24=dGf_Bi14OS24*F;
dGf_Bi2SO2=dGf_Bi2SO2*F;

dGf_AgBiS2=dGf_AgBiS2*F;
dGf_AgBi3S5=dGf_AgBi3S5*F;
dGf_AgBi6S9=dGf_AgBi6S9*F;
dGf_AgBiSO=dGf_AgBiSO*F;

dGf_H2O=  dGf_H2O*F;          
###############################################################################
###############################################################################
###############################################################################


###############################################################################
###############   Populate the species matrix  ################################
###############################################################################
species=np.zeros((59,8))

######## Formation Energies ###################################################
species[0,0]=0.00;
species[1,0]=dGf_AgO
species[2,0]=dGf_Ag2O
species[3,0]=dGf_Ag2O3
species[4,0]=dGf_Ag_Plus
species[5,0]=dGf_Ag_2Plus
species[6,0]=dGf_AgO_Plus
species[7,0]=dGf_AgO_Minus
species[8,0]=dGf_AgOH
species[9,0]=dGf_AgOH2_Minus

species[10,0]=dGf_Bi
species[11,0]=dGf_Bi2O3
species[12,0]=dGf_Bi2O5
species[13,0]=dGf_Bi2O4
species[14,0]=dGf_Bi4O7
species[15,0]=dGf_Bi_3Plus
species[16,0]=dGf_BiOH_2Plus
species[17,0]=dGf_BiO_Plus

species[18,0]=dGf_S
species[19,0]=dGf_H2S
species[20,0]=dGf_HS_Minus
species[21,0]=dGf_S_2Minus
species[22,0]=dGf_S2_2Minus
species[23,0]=dGf_S3_2Minus
species[24,0]=dGf_S4_2Minus
species[25,0]=dGf_S5_2Minus
species[26,0]=dGf_H2S2O3
species[27,0]=dGf_HS2O3_Minus
species[28,0]=dGf_S2O3_2Minus
species[29,0]=dGf_S5O6_2Minus
species[30,0]=dGf_S4O6_2Minus
species[31,0]=dGf_HS2O4_Minus
species[32,0]=dGf_S2O4_2Minus
species[33,0]=dGf_S3O6_2Minus
species[34,0]=dGf_H2SO3
species[35,0]=dGf_HSO3_Minus
species[36,0]=dGf_SO3_2Minus
species[37,0]=dGf_S2O6_2Minus
species[38,0]=dGf_H2SO4
species[39,0]=dGf_HSO4_Minus
species[40,0]=dGf_SO4_2Minus
species[41,0]=dGf_S2O8_2Minus
species[42,0]=dGf_HSO5_Minus
species[43,0]=dGf_S2O5_2Minus

species[44,0]=dGf_Ag2S
species[45,0]=dGf_Ag8S

species[46,0]=dGf_Ag5BiO4
species[47,0]=dGf_Ag2BiO3
species[48,0]=dGf_AgBiO3

species[49,0]=dGf_BiS2
species[50,0]=dGf_Bi2S3
species[51,0]=dGf_Bi2S2O
species[52,0]=dGf_Bi2SO4_3
species[53,0]=dGf_Bi14OS24
species[54,0]=dGf_Bi2SO2

species[55,0]=dGf_AgBiS2
species[56,0]=dGf_AgBi3S5
species[57,0]=dGf_AgBi6S9
species[58,0]=dGf_AgBiSO

######## Electron Count #######################################################
 #Ag
species[0,1]=0;
species[1,1]=2
species[2,1]=2
species[3,1]=6
species[4,1]=1
species[5,1]=2
species[6,1]=3
species[7,1]=1
species[8,1]=1
species[9,1]=1
#Bi
species[10,1]=0
species[11,1]=6
species[12,1]=10
species[13,1]=8
species[14,1]=14
species[15,1]=3
species[16,1]=3
species[17,1]=3
#S
species[18,1]=0
species[19,1]=-2
species[20,1]=-2
species[21,1]=-2
species[22,1]=-2
species[23,1]=-2
species[24,1]=-2
species[25,1]=-2
species[26,1]=4
species[27,1]=4
species[28,1]=4
species[29,1]=10
species[30,1]=10
species[31,1]=6
species[32,1]=6
species[33,1]=10
species[34,1]=4
species[35,1]=4
species[36,1]=4
species[37,1]=10
species[38,1]=6
species[39,1]=6
species[40,1]=6
species[41,1]=14  
species[42,1]=8
species[43,1]=8     
#AgBiOS
species[44,1]=0
species[45,1]=0
species[46,1]=8
species[47,1]=6
species[48,1]=6
#BiSO
species[49,1]=0
species[50,1]=0
species[51,1]=2
species[52,1]=24
species[53,1]=2
species[54,1]=4
#BiAgSO
species[55,1]=0
species[56,1]=0
species[57,1]=0
species[58,1]=2

######## Hydrogen H+ Count ####################################################
#Ag
species[0,2]=0
species[1,2]=2
species[2,2]=2
species[3,2]=6
species[4,2]=0
species[5,2]=0
species[6,2]=2
species[7,2]=2
species[8,2]=1
species[9,2]=2
#Bi
species[10,2]=0
species[11,2]=6
species[12,2]=10
species[13,2]=8
species[14,2]=14
species[15,2]=0
species[16,2]=1
species[17,2]=2
#S
species[18,2]=0
species[19,2]=-2
species[20,2]=-1
species[21,2]=0
species[22,2]=0
species[23,2]=0
species[24,2]=0
species[25,2]=0
species[26,2]=6 
species[27,2]=5
species[28,2]=4
species[29,2]=12
species[30,2]=12
species[31,2]=6
species[32,2]=8
species[33,2]=12
species[34,2]=4
species[35,2]=5
species[36,2]=6
species[37,2]=12
species[38,2]=6
species[39,2]=7
species[40,2]=8
species[41,2]=16
species[42,2]=9
species[43,2]=10
#AgSOBi
species[44,2]=0
species[45,2]=0
species[46,2]=8
species[47,2]=6
species[48,2]=6
#BiSO
species[49,2]=0
species[50,2]=0
species[51,2]=2
species[52,2]=24
species[53,2]=2
species[54,2]=4
#BiAgS
species[55,2]=0
species[56,2]=0
species[57,2]=0
species[58,2]=2


######## Number Silvers Ag  ###################################################
#Silver
species[0,3]=1
species[1,3]=1
species[2,3]=2
species[3,3]=2
species[4,3]=1
species[5,3]=1
species[6,3]=1
species[7,3]=1
species[8,3]=1
species[9,3]=1
#Bismuth
species[10,3]=0
species[11,3]=0
species[12,3]=0
species[13,3]=0
species[14,3]=0
species[15,3]=0
species[16,3]=0
species[17,3]=0
#Sulphur
species[18,3]=0
species[19,3]=0
species[20,3]=0
species[21,3]=0
species[22,3]=0
species[23,3]=0
species[24,3]=0
species[25,3]=0
species[26,3]=0
species[27,3]=0
species[28,3]=0
species[29,3]=0
species[30,3]=0
species[31,3]=0
species[32,3]=0
species[33,3]=0
species[34,3]=0
species[35,3]=0
species[36,3]=0
species[37,3]=0
species[38,3]=0
species[39,3]=0
species[40,3]=0
species[41,3]=0
species[42,3]=0
species[43,3]=0
#AgBiSO
species[44,3]=2
species[45,3]=8
species[46,3]=5
species[47,3]=2
species[48,3]=1
#BiSO
species[49,3]=0
species[50,3]=0
species[51,3]=0
species[52,3]=0
species[53,3]=0
species[54,3]=0
#AgBi
species[55,3]=1
species[56,3]=1
species[57,3]=1
species[58,3]=1

########### Number of Bismuths Bi #############################################
#Silver
species[0,4]=0
species[1,4]=0
species[2,4]=0
species[3,4]=0
species[4,4]=0
species[5,4]=0
species[6,4]=0
species[7,4]=0
species[8,4]=0
species[9,4]=0
#Bismuth
species[10,4]=1
species[11,4]=2
species[12,4]=2
species[13,4]=2
species[14,4]=4
species[15,4]=1
species[16,4]=1
species[17,4]=1
#Sulphur
species[18,4]=0
species[19,4]=0
species[20,4]=0
species[21,4]=0
species[22,4]=0
species[23,4]=0
species[24,4]=0
species[25,4]=0
species[26,4]=0
species[27,4]=0
species[28,4]=0
species[29,4]=0
species[30,4]=0
species[31,4]=0
species[32,4]=0
species[33,4]=0
species[34,4]=0
species[35,4]=0
species[36,4]=0
species[37,4]=0
species[38,4]=0
species[39,4]=0
species[40,4]=0
species[41,4]=0
species[42,4]=0
species[43,4]=0
#AgSBiO
species[44,4]=0
species[45,4]=0
species[46,4]=1
species[47,4]=1
species[48,4]=1
#BiSO
species[49,4]=1
species[50,4]=2
species[51,4]=2
species[52,4]=2
species[53,4]=14
species[54,4]=2
#BiAgSO
species[55,4]=1
species[56,4]=3
species[57,4]=6
species[58,4]=1

########### Number of Sulphurs S  #############################################
#Silver
species[0,5]=0
species[1,5]=0
species[2,5]=0
species[3,5]=0
species[4,5]=0
species[5,5]=0
species[6,5]=0
species[7,5]=0
species[8,5]=0
species[9,5]=0
#Bismuth
species[10,5]=0
species[11,5]=0
species[12,5]=0
species[13,5]=0
species[14,5]=0
species[15,5]=0
species[16,5]=0
species[17,5]=0
#Sulphur
species[18,5]=1
species[19,5]=1
species[20,5]=1
species[21,5]=1
species[22,5]=2
species[23,5]=3
species[24,5]=4
species[25,5]=5
species[26,5]=2
species[27,5]=2
species[28,5]=2
species[29,5]=5
species[30,5]=4
species[31,5]=2
species[32,5]=2
species[33,5]=3
species[34,5]=1
species[35,5]=1
species[36,5]=1
species[37,5]=2
species[38,5]=1
species[39,5]=1
species[40,5]=1
species[41,5]=2
species[42,5]=1
species[43,5]=2
#AgBiSO
species[44,5]=1
species[45,5]=1
species[46,5]=0
species[47,5]=0
species[48,5]=0
#BiSO
species[49,5]=2
species[50,5]=3
species[51,5]=2
species[52,5]=3
species[53,5]=24
species[54,5]=1
#BiAgSO
species[55,5]=2
species[56,5]=5
species[57,5]=9
species[58,5]=1

######### Number of H2O's #####################################################
#Silver
species[0,6]=0
species[1,6]=1
species[2,6]=1
species[3,6]=3
species[4,6]=0
species[5,6]=0
species[6,6]=1
species[7,6]=1
species[8,6]=1
species[9,6]=2
#Bi
species[10,6]=0
species[11,6]=3
species[12,6]=5
species[13,6]=4
species[14,6]=7
species[15,6]=0
species[16,6]=1
species[17,6]=1
#Sulphur
species[18,6]=0
species[19,6]=0
species[20,6]=0
species[21,6]=0
species[22,6]=0
species[23,6]=0
species[24,6]=0
species[25,6]=0
species[26,6]=3
species[27,6]=3
species[28,6]=3
species[29,6]=6
species[30,6]=6
species[31,6]=4
species[32,6]=4
species[33,6]=6
species[34,6]=3
species[35,6]=3
species[36,6]=3
species[37,6]=6
species[38,6]=4
species[39,6]=4
species[40,6]=4
species[41,6]=8
species[42,6]=5
species[43,6]=5
#AgBiSO
species[44,6]=0
species[45,6]=0
species[46,6]=4
species[47,6]=3
species[48,6]=3
#BiSO
species[49,6]=0
species[50,6]=0
species[51,6]=1
species[52,6]=12
species[53,6]=1
species[54,6]=2
#AgBiSO
species[55,6]=0
species[56,6]=0
species[57,6]=0
species[58,6]=1


########## Aqueous Ions?????? #################################################
#Silver
species[0,7]=0
species[1,7]=0
species[2,7]=0
species[3,7]=0
species[4,7]=1
species[5,7]=1
species[6,7]=1
species[7,7]=1
species[8,7]=1
species[9,7]=1
#Bismuth
species[10,7]=0
species[11,7]=0
species[12,7]=0
species[13,7]=0
species[14,7]=0
species[15,7]=1
species[16,7]=1
species[17,7]=1
#Sulphur
species[18,7]=0
species[19,7]=1
species[20,7]=1
species[21,7]=1
species[22,7]=1
species[23,7]=1
species[24,7]=1
species[25,7]=1
species[26,7]=1
species[27,7]=1
species[28,7]=1
species[29,7]=1
species[30,7]=1
species[31,7]=1
species[32,7]=1
species[33,7]=1
species[34,7]=1
species[35,7]=1
species[36,7]=1
species[37,7]=1
species[38,7]=1
species[39,7]=1
species[40,7]=1
species[41,7]=1
species[42,7]=1
species[43,7]=1
#AgSBiO
species[44,7]=0
species[45,7]=0
species[46,7]=0
species[47,7]=0
species[48,7]=0
#BiSO
species[49,7]=0
species[50,7]=0
species[51,7]=0
species[52,7]=0
species[53,7]=0
species[54,7]=0
#AgBiSO
species[55,7]=0
species[56,7]=0
species[57,7]=0
species[58,7]=0

#Function to determine species combinations
try:
    combos=load('BiAgOS-speciesCombo.npy')
    num=load('BiAgOS-numberSpecies.npy')
    combo_num=int(num[0])
except OSError:
    print('Cannot Open File')
    ###############################################################################
    #### Determine which species are able to combine at the composition ###########
    ###############################################################################
    t=1
    flag=1
    f_total=int;
    f=np.zeros((3))
    combos=np.zeros((40000,9,3))
    combo_num=0
    combos[combo_num, 0, 0]=-1
    combos[combo_num, 0, 1]=-1
    combos[combo_num, 0, 2]=-1
    for k in range(0, len(species)):
        for m in range(0, len(species)):
            for p in range(0, len(species)):
                #Check to make sure each element is in this combination of species
                if((species[k, 3]>0 or species[m, 3] >0 or species[p, 3]) \
                   and (species[k, 4]>0 or species[m, 4] >0 or species[p, 4]) \
                   and (species[k, 5]>0 or species[m, 5] >0 or species[p, 5])):
                    #save species in array
                    t=1
                    a = np.array([[species[k, 3],species[m, 3], species[p,3]], \
                                  [species[k, 4],species[m, 4], species[p,4]], \
                                  [species[k, 5],species[m, 5], species[p,5]]])
                        
                    #check to see if each species contains a single element. This is a really long call.
                    flag=1
                    if((species[k, 3]==0 and species[m, 3] ==0) or \
                       (species[m, 3]==0 and species[p, 3] ==0) or \
                       (species[k, 3]==0 and species[p, 3] ==0)):
                        if((species[k, 4]==0 and species[m, 4] ==0) or \
                           (species[m, 4]==0 and species[p, 4] ==0) or \
                           (species[k, 4]==0 and species[p, 4] ==0)):
                            if((species[k, 5]==0 and species[m, 5] ==0) or \
                               (species[m, 5]==0 and species[p, 5] ==0) or \
                               (species[k, 5]==0 and species[p, 5] ==0)):
                
                                flag=0
                                #if so, find the composition though linear algebra.
                                try:
                                    f=np.linalg.solve(a, composition)
                                except:
                                    #print('Error: Species '+str(k)+', Species2: '+str(m)+', Species3: '+str(p)+'\n') 
                                    t=1
                                t=0
                    

                    #If there is at least one multi-element species in this combination
                    if(flag==1):
                        #test each linear combination
                        for h in range(1, 20):
                            for i in range(1, 20):
                                for j in range(1, 20):
                                    #Is there a linear combination of the elements that will allow 
                                    #For the 
                                    if(((h*a[0,0]+i*a[0,1]+j*a[0,2])/(h*a[1,0]+i*a[1,1]+j*a[1,2]))==composition[0]/composition[1] and \
                                       ((h*a[1,0]+i*a[1,1]+j*a[1,2])/(h*a[2,0]+i*a[2,1]+j*a[2,2]))==composition[1]/composition[2] and \
                                       ((h*a[0,0]+i*a[0,1]+j*a[0,2])/(h*a[2,0]+i*a[2,1]+j*a[2,2]))==composition[0]/composition[2]):
                                        #save the composition
                                        f[0]=h
                                        f[1]=i
                                        f[2]=j
                                        #Ending parameters, break loops
                                        t=0;
                                        h=40;
                                        i=40;
                                        j=40;
                            
                                    #If there is a linear combination, save the species in the combos array.
            
                    if (t==0):
                        #print(str(combo_num)+': Species1: '+str(k)+', Species2: '+str(m)+'\n')
                        #Species Number
                        combos[combo_num, 0, 0]=k
                        combos[combo_num, 0, 1]=m
                        combos[combo_num, 0, 2]=p
                        #Energy
                        combos[combo_num, 1, 0]=species[k,0]
                        combos[combo_num, 1, 1]=species[m,0]
                        combos[combo_num, 1, 2]=species[p,0]
                        #Electrons
                        combos[combo_num, 2, 0]=species[k,1]
                        combos[combo_num, 2, 1]=species[m,1]
                        combos[combo_num, 2, 2]=species[p,1]
                        #H+
                        combos[combo_num, 3, 0]=species[k,2]
                        combos[combo_num, 3, 1]=species[m,2]
                        combos[combo_num, 3, 2]=species[p,2]
                        #Number Silvers
                        combos[combo_num, 4, 0]=species[k,3]
                        combos[combo_num, 4, 1]=species[m,3]
                        combos[combo_num, 4, 2]=species[p,3]
                        #Number Bismuth
                        combos[combo_num, 5, 0]=species[k,4]
                        combos[combo_num, 5, 1]=species[m,4]
                        combos[combo_num, 5, 2]=species[p,4]
                        #Number H2O
                        combos[combo_num, 6, 0]=species[k,5]
                        combos[combo_num, 6, 1]=species[m,5]
                        combos[combo_num, 6, 2]=species[p,5]
                        #Aqueous Ions
                        combos[combo_num, 7, 0]=species[k,6]
                        combos[combo_num, 7, 1]=species[m,6]
                        combos[combo_num, 7, 2]=species[p,6]
                        #Percent of each in species in final combo
                        f_total=f[0]+f[1]+f[2];
                        combos[combo_num, 8, 0]=f[0]/f_total
                        combos[combo_num, 8, 1]=f[1]/f_total
                        combos[combo_num, 8, 2]=f[2]/f_total
                        
                        combo_num=combo_num+1;
                        t=1
                        #print('entered')
                    else:
                        #Catch and switch the value of t back to no
                        t=1

    save('BiAgOS-speciesCombo.npy', combos)
    save('BiAgOS-numberSpecies.npy', asarray([[combo_num]]))
    print('The number of species combinations is '+ str(combo_num)+'.\n')
###############################################################################
###############################################################################
###############################################################################


###############################################################################
###########  Chemical Potential Mesh Calculations  ############################
###############################################################################
#should be as long as there are specicies considered
#populate with smaller values that will be calculated.
muValues=np.zeros((n+1,n+1,4))
current_mu=int
current_ele=int
current_H=int
current_H2O=int
current_aquI=int
current_NumEle=int
sort=np.zeros((3,1))

#fill in the grid. Calculate 
for i in range(0, n+1):
    #calculate the energies for each species number
    pH=lowpH+(i*pHcount);
    for j in range(0,n+1):
        U=Ulow+(j*Ucount);
        muValues[i,j,0]=-1
        muValues[i,j,1]=-1
        muValues[i,j,2]=-1
        muValues[i,j,3]=100000000
        #Go through all species, commpare all pairs
        for k in range(0, combo_num):
            p=int(combos[k,0,0]);
            m=int(combos[k,0,1]);
            s=int(combos[k,0,2]);
            f1=combos[k,8,0];
            f2=combos[k,8,1];
            f3=combos[k,8,2];
            
            #The first species's contribution to the mu
            current_eng=species[p,0]
            current_ele=F*U*(species[p,1])
            current_H=R*T*np.log(10.0)*pH*(species[p,2])
            current_H2O=dGf_H2O*(species[p,6])
            current_aquI=R*T*np.log(nI)*(species[p,7])
            current_NumEle=1
            for t in range(3,6):
                if(species[p,t]>1):
                    current_NumEle=current_NumEle*species[p,t];
            current_mu=f1*((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
             
            #The second species' contribution to the mu
            current_eng=species[m,0];
            current_ele=F*U*(species[m,1])
            current_H=R*T*np.log(10.0)*pH*(species[m,2])
            current_H2O=dGf_H2O*(species[m,6])
            current_aquI=R*T*np.log(nI)*(species[m,7])
            current_NumEle=1
            for t in range(3,6):
                if(species[m,t]>1):
                    current_NumEle=current_NumEle*species[m,t];
            current_mu=current_mu+f2*((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
            
            #The second species' contribution to the mu
            current_eng=species[s,0];
            current_ele=F*U*(species[s,1])
            current_H=R*T*np.log(10.0)*pH*(species[s,2])
            current_H2O=dGf_H2O*(species[s,6])
            current_aquI=R*T*np.log(nI)*(species[s,7])
            current_NumEle=1
            for t in range(3,6):
                if(species[s,t]>1):
                    current_NumEle=current_NumEle*species[s,t];
            current_mu=current_mu+f3*((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
            if(current_mu<muValues[i,j,3]):
                sort[0,0]=p
                sort[1,0]=m
                sort[2,0]=s
                a=np.sort(sort[:,0])
                muValues[i,j,0]=a[0]
                muValues[i,j,1]=a[1]
                muValues[i,j,2]=a[2]
                muValues[i,j,3]=current_mu
###############################################################################
###############################################################################
###############################################################################
                
                
            
###############################################################################
###################  Plot Pourbaix Diagram  ###################################
###############################################################################
flag = np.zeros((50,6)) # The first 4 indexes are the materials stored, the next three are the colors
index=0;
fig =plt.figure()
ax=plt.subplot(111)
ax = plt.gca() 
ax.set_xlim([lowpH,highpH])
ax.set_ylim([Ulow,Uhigh])     
l=0;
index=0;
for i in range(0, n+1):
    pH=lowpH+i*pHcount;
    for j in range(0,n+1):
        U=Ulow+(Ucount*j);
        l=0
        for k in range(0, len(flag)):
            if(flag[k,0]==muValues[i,j,0] and flag[k,1]==muValues[i,j,1] and flag[k,2]==muValues[i,j,2]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,0] and flag[k,1]==muValues[i,j,2]and flag[k,2]==muValues[i,j,1]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,1] and flag[k,1]==muValues[i,j,2]and flag[k,2]==muValues[i,j,0]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,1] and flag[k,1]==muValues[i,j,0]and flag[k,2]==muValues[i,j,2]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,2] and flag[k,1]==muValues[i,j,0]and flag[k,2]==muValues[i,j,1]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,2] and flag[k,1]==muValues[i,j,1]and flag[k,2]==muValues[i,j,0]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
        if(l==0):
            label='M1: '+str(muValues[i,j,0])+', M2: '+str(muValues[i,j,1])+' M3: '+str(muValues[i,j,2])
            flag[index,0] = muValues[i,j,0]
            flag[index,1] = muValues[i,j,1]
            flag[index,2] = muValues[i,j,2]
            flag[index,3] = random.random();
            flag[index,4] = random.random();
            flag[index,5] = random.random();
            ax.plot(pH,U,'.', color = [flag[index,3],flag[index,4],flag[index,5]],markersize=4,label=label)
            index=index+1;


#####Plot H2O and H2 lines##################################
muH=np.zeros((pHrange+1));
muH2O=np.zeros((pHrange+1));
pHArray=np.zeros((pHrange+1));
for i in range(0, pHrange):
    pHArray[i] =lowpH+i;
    muH[i]=-0.059*pHArray[i];
    muH2O[i]=1.23-0.059*pHArray[i];

pHArray[pHrange] =lowpH+(pHrange);
muH[pHrange]=-0.059*pHArray[pHrange];
muH2O[pHrange]=1.23-0.059*pHArray[pHrange];
##############################################################
            
ax.plot(pHArray[:], muH[:],'c--',label='$H_2$',linewidth=1)
ax.plot(pHArray[:], muH2O[:],'b--',label='$H_2O$', linewidth=1)            
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 0.9), ncol=1)
plt.ylabel('Electric Potential, E(V)')
plt.xlabel('pH')
plt.title('Bi-Ag-S Pourbaix Diagram, $\eta_{Bi,Ag,S}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Ag:' +str(composition[1])+'Bi:'+str(composition[2])+'S')


###############################################################################
##############   Plot with Lines   ############################################
###############################################################################
flag = np.zeros((50,6)) # The first 4 indexes are the materials stored, the next three are the colors
index=0;
fig =plt.figure()
ax=plt.subplot(111)
ax = plt.gca() 
ax.set_xlim([lowpH,highpH])
ax.set_ylim([Ulow,Uhigh]) 

#If drawing lines for metastable phases
for i in range(1, n):
    #calculate the energies for each species number
    pH=lowpH+(i*pHcount);
    for j in range(1,n):
        U=Ulow+(j*Ucount);      
        #If drawing lines for metastable phases
        
        if((muValues[i,j,0]!=muValues[i-1,j,0])):
            ax.plot(pH,U,'.', color = [0.0,0.0,0.0],markersize=2)
        elif(muValues[i,j,1]!=muValues[i-1,j,1]):
            ax.plot(pH,U,'.', color = [0.0,0.0,0.0],markersize=2)
        elif((muValues[i,j,0]!=muValues[i,j-1,0]) or (muValues[i,j,1]!=muValues[i,j-1,1])):
            ax.plot(pH,U,'.', color = [0.0,0.0,0.0],markersize=2)
        elif((muValues[i,j,2]!=muValues[i,j-1,2]) or (muValues[i,j,2]!=muValues[i-1,j,2])):
            ax.plot(pH,U,'.', color = [0.0,0.0,0.0],markersize=2)

ax.plot(pHArray[:], muH[:],'c--',label='$H_2$',linewidth=1)
ax.plot(pHArray[:], muH2O[:],'b--',label='$H_2O$', linewidth=1)  
        
plt.ylabel('Electric Potential, E(V)')
plt.xlabel('pH')
plt.title('Bi-Ag-S Pourbaix Diagram, $\eta_{Bi,Ag,S}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Ag:' +str(composition[1])+'Bi:'+str(composition[2])+'S')
chartBox=ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.5, chartBox.height*1.5])
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 0.9), ncol=1)
plt.show()



print('End of Script')