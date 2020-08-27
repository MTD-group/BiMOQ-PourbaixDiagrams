#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:17:41 2018

@author: laurenwalters
"""
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

#Array showing the composition of Cu:Bi:S
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
Ee_Cu=-4.3152965;        
Ee_CuO=-10.7488868;
Ee_Cu2O=-14.99698;
Ee_CuOH2_s=-25.1916025;
#PBEsol
Ee_O2=-10.281123
Ee_H2=-6.5141508
Ee_S= -4.391811875;
###############################################################################
########### MULTICOMPONENT SPECIES ############################################
#Calculated with PBEsol
Ee_Cu2S=-13.4793116667;
Ee_Cu7S4=-49.8241325;
Ee_CuS=-9.170266;
Ee_CuS2=-13.63935;
Ee_Cu2SO4_3=-101.5166;

Ee_BiCu=-9.31218;
Ee_CuBiO2_2=-42.245475;

Ee_BiS2=-14.6172585;
Ee_Bi2S3=-24.878388;
Ee_Bi2S2O=-27.2327565;
Ee_Bi2SO4_3=-109.35902;
Ee_Bi14OS24=-247.57619;
Ee_Bi2SO2=-29.50652;

Ee_BiSCuO=-21.5022935;
Ee_Cu3BiS3=-32.4713275;
Ee_Cu4Bi4S9=-80.830705;
Ee_Cu4BiS2_5=-90.647798;
Ee_CuBiS2=-19.041996;
###############################################################################
###### Vibrational Energy #####################################################
###############################################################################
#Vibrational Energies in eV/f.u.
#From PBEsol Phonon Calculations
Fvib_O2=-0.272;
F_rot_trans_O2=0.099;
Ftot_O2=Fvib_O2+F_rot_trans_O2;
F_H = .202;

Fvib_S=-0.0091266451372

Fvib_CuO=0.062498987735
Fvib_Cu2O=0.00507624852
Fvib_Cu=-0.007167374680
Fvib_CuOH2_s=0.66653026525

Fvib_Bi=-0.0761976993239
Fvib_Bi2O3=-0.057653546889
Fvib_Bi2O5=0.14677315404
Fvib_Bi2O4=0.12231438709
Fvib_Bi4O7=0.08741679245

Fvib_Cu2S=-0.0050937891364
Fvib_Cu7S4=-0.178002185722
Fvib_CuS=-0.0119849701814
Fvib_CuS2=-0.0033060080158
Fvib_Cu2SO4_3=1.00135494361

Fvib_BiCu=-0.11006963132
Fvib_CuBiO2_2=0.09853363658

Fvib_BiS2=-0.063943629448
Fvib_Bi2S3=-0.1428187610337
Fvib_Bi2S2O=-0.08193190191
Fvib_Bi2SO4_3=0.81266278392
Fvib_Bi14OS24=0.02990373431
Fvib_Bi2SO2=-0.0265520338422

Fvib_BiSCuO=-0.039894146059
Fvib_Cu3BiS3=-0.1661179102334
Fvib_Cu4Bi4S9=-0.3270592722135
Fvib_Cu4BiS2_5=-0.430548296696
Fvib_CuBiS2=-0.08663072302
###############################################################################
###   Compounds-Calculate the formation energies   ############################
###############################################################################
#Free Energies of Formation in eV/f.u.
dGf_CuO= (Ee_CuO+Fvib_CuO)             -(Ee_Cu+Fvib_Cu) - 0.5*(Ee_O2+Ftot_O2);
dGf_Cu2O=(Ee_Cu2O+Fvib_Cu2O)           -2.0*(Ee_Cu+Fvib_Cu) - 0.5*(Ee_O2+Ftot_O2);
dGf_CuOH2_s= (Ee_CuOH2_s+Fvib_CuOH2_s) -(Ee_Cu+Fvib_Cu)-(Ee_O2+Ftot_O2)-(Ee_H2+F_H);

dGf_Bi2O3= ((Ee_Bi2O3)+Fvib_Bi2O3) -2.0*(Ee_Bi+Fvib_Bi)-1.5*(Ee_O2-Ftot_O2); 
dGf_Bi2O5= ((Ee_Bi2O5)+Fvib_Bi2O5) -2.0*(Ee_Bi+Fvib_Bi)-2.5*(Ee_O2-Ftot_O2);
dGf_Bi2O4= ((Ee_Bi2O4)+Fvib_Bi2O4) -2.0*(Ee_Bi+Fvib_Bi)-2.0*(Ee_O2-Ftot_O2);
dGf_Bi4O7= ((Ee_Bi4O7)+Fvib_Bi4O7) -4.0*(Ee_Bi+Fvib_Bi)-3.5*(Ee_O2-Ftot_O2);

dGf_Cu2S=(Ee_Cu2S+Fvib_Cu2S)             -2*(Ee_Cu+Fvib_Cu)-(Ee_S+Fvib_S);
dGf_Cu7S4=(Ee_Cu7S4+Fvib_Cu7S4)          -7*(Ee_Cu+Fvib_Cu)-4*(Ee_S+Fvib_S);
dGf_CuS=(Ee_CuS+Fvib_CuS)                -(Ee_Cu+Fvib_Cu)-(Ee_S+Fvib_S);
dGf_CuS2=(Ee_CuS2+Fvib_CuS2)             -(Ee_Cu+Fvib_Cu)-2*(Ee_S+Fvib_S);
dGf_Cu2SO4_3=(Ee_Cu2SO4_3+Fvib_Cu2SO4_3) -2*(Ee_Cu+Fvib_Cu)-3*(Ee_S+Fvib_S)-6.0*((Ee_O2)-Ftot_O2);


dGf_BiCu=(Ee_BiCu+Fvib_BiCu)              -(Ee_Cu+Fvib_Cu)-(Ee_Bi+Fvib_Bi);
dGf_CuBiO2_2=(Ee_CuBiO2_2+Fvib_CuBiO2_2)  -(Ee_Cu+Fvib_Cu)-2*(Ee_Bi+Fvib_Bi)-2.0*((Ee_O2)-Ftot_O2);

dGf_BiS2=(Ee_BiS2+Fvib_BiS2)                  -(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S);
dGf_Bi2S3=(Ee_Bi2S3+Fvib_Bi2S3)               -2*(Ee_Bi+Fvib_Bi)-3*(Ee_S+Fvib_S);
dGf_Bi2S2O=(Ee_Bi2S2O+Fvib_Bi2S2O)            -2*(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);
dGf_Bi2SO4_3=(Ee_Bi2SO4_3+Fvib_Bi2SO4_3)      -2*(Ee_Bi+Fvib_Bi)-3*(Ee_S+Fvib_S)-6.0*((Ee_O2)-Ftot_O2);
dGf_Bi14OS24=(Ee_Bi14OS24+Fvib_Bi14OS24)      -14*(Ee_Bi+Fvib_Bi)-24*(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);
dGf_Bi2SO2=(Ee_Bi2SO2+Fvib_Bi2SO2)            -2*(Ee_Bi+Fvib_Bi)-(Ee_S+Fvib_S)-1.0*((Ee_O2)-Ftot_O2);

dGf_BiSCuO=(Ee_BiSCuO+Fvib_BiSCuO)         -(Ee_Cu+Fvib_Cu)-(Ee_Bi+Fvib_Bi)-(Ee_S+Fvib_S)-0.5*((Ee_O2)-Ftot_O2);
dGf_Cu3BiS3=(Ee_Cu3BiS3+Fvib_Cu3BiS3)      -3*(Ee_Cu+Fvib_Cu)-(Ee_Bi+Fvib_Bi)-3*(Ee_S+Fvib_S);
dGf_Cu4Bi4S9=(Ee_Cu4Bi4S9+Fvib_Cu4Bi4S9)   -4*(Ee_Cu+Fvib_Cu)-4*(Ee_Bi+Fvib_Bi)-9*(Ee_S+Fvib_S);
dGf_Cu4BiS2_5=(Ee_Cu4BiS2_5+Fvib_Cu4BiS2_5)-4*(Ee_Cu+Fvib_Cu)-5*(Ee_Bi+Fvib_Bi)-10*(Ee_S+Fvib_S);
dGf_CuBiS2=(Ee_CuBiS2+Fvib_CuBiS2)         -(Ee_Cu+Fvib_Cu)-(Ee_Bi+Fvib_Bi)-2*(Ee_S+Fvib_S);

#Set the reference values
dGf_Cu=0.0;
dGf_Bi=0.0;
dGf_S=0.0;
###############################################################################
###############################################################################
###############################################################################


###############################################################################
############## Aqueous Ion Free Energies of Formation #########################
#Free Energies of Formation in eV/f.u.
##Elemental Bismuth Species
dGf_Bi_3Plus= 0.6430898
dGf_BiOH_2Plus=  -1.6968378
dGf_BiO_Plus=  -1.4977965

##Elemental Copper Species
dGf_Cu1=  0.506502
dGf_Cu2=  0.674092
dGf_CuOH2_minus= -3.4518209
dGf_CuOH3= -5.1197432
dGf_CuOH_Plus= -1.3127387
dGf_CuOH4_2=-6.814302
dGf_CuOH2= -3.2666113 
dGf_CuOH = -1.2677578
dGf_Cu2OH2_2plus=-2.942417
dGf_Cu3OH4_2plus=-6.567839

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

################################################################################
#############   CONVERT from eV to kJ/mol  ####################################
###############################################################################
dGf_Cu=   dGf_Cu*F;
dGf_CuO=  dGf_CuO*F;
dGf_Cu2O= dGf_Cu2O*F;
dGf_Cu1=  dGf_Cu1*F;
dGf_Cu2=  dGf_Cu2*F;
dGf_CuOH4_2= dGf_CuOH4_2*F;
dGf_CuOH2_minus= dGf_CuOH2_minus*F;
dGf_CuOH3= dGf_CuOH3*F;
dGf_CuOH_Plus= dGf_CuOH_Plus*F;
dGf_CuOH2= dGf_CuOH2*F;
dGf_CuOH = dGf_CuOH*F;
dGf_Cu2OH2_2plus=dGf_Cu2OH2_2plus*F;
dGf_Cu3OH4_2plus=dGf_Cu3OH4_2plus*F;
dGf_CuOH2_s=dGf_CuOH2_s*F

dGf_Bi=   dGf_Bi*F; 
dGf_Bi2O3=  dGf_Bi2O3*F; 
dGf_Bi2O5= dGf_Bi2O5*F; 
dGf_Bi2O4=dGf_Bi2O4*F;
dGf_Bi4O7=dGf_Bi4O7*F;
dGf_Bi_3Plus= dGf_Bi_3Plus*F;
dGf_BiOH_2Plus=  dGf_BiOH_2Plus*F;
dGf_BiO_Plus=  dGf_BiO_Plus*F;

dGf_S=   dGf_S*F;              
dGf_H2S=dGf_H2S*F;             
dGf_HS_Minus=dGf_HS_Minus*F;         
dGf_S_2Minus=dGf_S_2Minus*F;   
dGf_S2_2Minus=dGf_S2_2Minus*F;    
dGf_S3_2Minus=dGf_S3_2Minus*F;
dGf_S4_2Minus=dGf_S4_2Minus*F;       
dGf_S5_2Minus=dGf_S5_2Minus*F;       
dGf_H2S2O3=dGf_H2S2O3*F;       
dGf_HS2O3_Minus=dGf_HS2O3_Minus*F;    
dGf_S2O3_2Minus=dGf_S2O3_2Minus*F;      
dGf_S5O6_2Minus=dGf_S5O6_2Minus*F;      
dGf_S4O6_2Minus=dGf_S4O6_2Minus*F;
dGf_HS2O4_Minus=dGf_HS2O4_Minus*F;  
dGf_S2O4_2Minus=dGf_S2O4_2Minus*F;    
dGf_S3O6_2Minus=dGf_S3O6_2Minus*F; 
dGf_H2SO3=dGf_H2SO3*F;         
dGf_HSO3_Minus=dGf_HSO3_Minus*F;      
dGf_SO3_2Minus=dGf_SO3_2Minus*F;       
dGf_S2O6_2Minus=dGf_S2O6_2Minus*F;        
dGf_H2SO4=dGf_H2SO4*F;         
dGf_HSO4_Minus=dGf_HSO4_Minus*F;  
dGf_SO4_2Minus=dGf_SO4_2Minus*F;  
dGf_S2O8_2Minus=dGf_S2O8_2Minus*F;       
dGf_HSO5_Minus=dGf_HSO5_Minus*F;
dGf_S2O5_2Minus=dGf_S2O5_2Minus*F;

dGf_Cu2S=dGf_Cu2S*F;
dGf_Cu7S4=dGf_Cu7S4*F;
dGf_CuS=dGf_CuS*F;
dGf_CuS2=dGf_CuS2*F;
dGf_Cu2SO4_3=dGf_Cu2SO4_3*F;


dGf_BiCu=dGf_BiCu*F;
dGf_CuBiO2_2=dGf_CuBiO2_2*F;

dGf_BiS2=dGf_BiS2*F;
dGf_Bi2S3=dGf_Bi2S3*F;
dGf_Bi2S2O=dGf_Bi2S2O*F;
dGf_Bi2SO4_3=dGf_Bi2SO4_3*F;
dGf_Bi14OS24=dGf_Bi14OS24*F;
dGf_Bi2SO2=dGf_Bi2SO2*F;

dGf_BiSCuO=dGf_BiSCuO*F;
dGf_Cu3BiS3=dGf_Cu3BiS3*F;
dGf_Cu4Bi4S9=dGf_Cu4Bi4S9*F;
dGf_Cu4BiS2_5=dGf_Cu4BiS2_5*F;
dGf_CuBiS2=dGf_CuBiS2*F;

dGf_H2O=  dGf_H2O*F;           
###############################################################################
###############################################################################
###############################################################################

###############################################################################
###############   Populate the species matrix  ################################
###############################################################################
species=np.zeros((65,8))

######## Formation Energies ###################################################
species[0,0]=0.00;
species[1,0]=dGf_CuO
species[2,0]=dGf_Cu2O
species[3,0]=dGf_Cu1
species[4,0]=dGf_Cu2
species[5,0]=dGf_CuOH4_2
species[6,0]=dGf_CuOH2_minus
species[7,0]=dGf_CuOH3
species[8,0]=dGf_CuOH_Plus
species[9,0]=dGf_CuOH2
species[10,0]=dGf_CuOH
species[11,0]=dGf_Cu2OH2_2plus
species[12,0]=dGf_Cu3OH4_2plus

species[13,0]=dGf_Bi
species[14,0]=dGf_Bi2O3
species[15,0]=dGf_Bi2O5
species[16,0]=dGf_Bi2O4
species[17,0]=dGf_Bi4O7
species[18,0]=dGf_Bi_3Plus
species[19,0]=dGf_BiOH_2Plus
species[20,0]=dGf_BiO_Plus

species[21,0]=dGf_S
species[22,0]=dGf_H2S
species[23,0]=dGf_HS_Minus
species[24,0]=dGf_S_2Minus
species[25,0]=dGf_S2_2Minus
species[26,0]=dGf_S3_2Minus
species[27,0]=dGf_S4_2Minus
species[28,0]=dGf_S5_2Minus
species[29,0]=dGf_H2S2O3
species[30,0]=dGf_HS2O3_Minus
species[31,0]=dGf_S2O3_2Minus
species[32,0]=dGf_S5O6_2Minus
species[33,0]=dGf_S4O6_2Minus
species[34,0]=dGf_HS2O4_Minus
species[35,0]=dGf_S2O4_2Minus
species[36,0]=dGf_S3O6_2Minus
species[37,0]=dGf_H2SO3
species[38,0]=dGf_HSO3_Minus
species[39,0]=dGf_SO3_2Minus
species[40,0]=dGf_S2O6_2Minus
species[41,0]=dGf_H2SO4
species[42,0]=dGf_HSO4_Minus
species[43,0]=dGf_SO4_2Minus
species[44,0]=dGf_S2O8_2Minus
species[45,0]=dGf_HSO5_Minus
species[46,0]=dGf_S2O5_2Minus

species[47,0]=dGf_Cu2S
species[48,0]=dGf_Cu7S4
species[49,0]=dGf_CuS
species[50,0]=dGf_CuS2
species[51,0]=dGf_Cu2SO4_3
species[52,0]=dGf_BiCu
species[53,0]=dGf_CuBiO2_2
species[54,0]=dGf_BiS2
species[55,0]=dGf_Bi2S3
species[56,0]=dGf_Bi2S2O
species[57,0]=dGf_Bi2SO4_3
species[58,0]=dGf_Bi14OS24
species[59,0]=dGf_Bi2SO2

species[60,0]=dGf_CuBiS2
species[61,0]=dGf_Cu4Bi4S9
species[62,0]=dGf_Cu4BiS2_5
species[63,0]=dGf_BiSCuO
species[64,0]=dGf_Cu3BiS3

######## Electron Count #######################################################
#Cu
species[0,1]=0.00;
species[1,1]=2
species[2,1]=2
species[3,1]=1
species[4,1]=2
species[5,1]=2
species[6,1]=1
species[7,1]=2
species[8,1]=2
species[9,1]=2
species[10,1]=1
species[11,1]=4
species[12,1]=6
#Bi
species[13,1]=0
species[14,1]=6
species[15,1]=10
species[16,1]=8
species[17,1]=14
species[18,1]=3
species[19,1]=3
species[20,1]=3
#S
species[21,1]=0
species[22,1]=-2
species[23,1]=-2
species[24,1]=-2
species[25,1]=-2
species[26,1]=-2
species[27,1]=-2
species[28,1]=-2
species[29,1]=4
species[30,1]=4
species[31,1]=4
species[32,1]=10
species[33,1]=10
species[34,1]=6
species[35,1]=6
species[36,1]=10
species[37,1]=4
species[38,1]=4
species[39,1]=4
species[40,1]=10
species[41,1]=6
species[42,1]=6
species[43,1]=6
species[44,1]=14
species[45,1]=8
species[46,1]=8        
#CuSOBi
species[47,1]=0
species[48,1]=0
species[49,1]=0
species[50,1]=0
species[51,1]=24
species[52,1]=0
species[53,1]=8
#BiSO
species[54,1]=0
species[55,1]=0
species[56,1]=2
species[57,1]=24
species[58,1]=2
species[59,1]=4
#CuBiS
species[60,1]=0
species[61,1]=0
species[62,1]=0
#BiCuSO
species[63,1]=2
species[64,1]=0

######## Hydrogen H+ Count ####################################################
#Cu
species[0,2]=0
species[1,2]=2
species[2,2]=2
species[3,2]=0
species[4,2]=0
species[5,2]=4
species[6,2]=2
species[7,2]=3
species[8,2]=1
species[9,2]=2
species[10,2]=1
species[11,2]=2
species[12,2]=4
#Bi
species[13,2]=0
species[14,2]=6
species[15,2]=10
species[16,2]=8
species[17,2]=14
species[18,2]=0
species[19,2]=1
species[20,2]=2
#S
species[21,2]=0
species[22,2]=-2
species[23,2]=-1
species[24,2]=0
species[25,2]=0
species[26,2]=0
species[27,2]=0
species[28,2]=0
species[29,2]=6
species[30,2]=5
species[31,2]=4
species[32,2]=12
species[33,2]=12
species[34,2]=6
species[35,2]=8
species[36,2]=12
species[37,2]=4
species[38,2]=5
species[39,2]=6
species[40,2]=12
species[41,2]=6
species[42,2]=7
species[43,2]=8
species[44,2]=16
species[45,2]=9
species[46,2]=10
#CuSBiO
species[47,2]=0
species[48,2]=0
species[49,2]=0
species[50,2]=0
species[51,2]=24
species[52,2]=0
species[53,2]=8     
#BiSO
species[54,2]=0
species[55,2]=0
species[56,2]=2
species[57,2]=24
species[58,2]=2
species[59,2]=4
#BiCuS
species[60,2]=0
species[61,2]=0
species[62,2]=0
#BiCuSO
species[63,2]=2
species[64,2]=0

########### Number of Coppers Cu ##############################################
#Cu
species[0,3]=1
species[1,3]=1
species[2,3]=2
species[3,3]=1
species[4,3]=1
species[5,3]=1
species[6,3]=1
species[7,3]=1
species[8,3]=1
species[9,3]=1
species[10,3]=1
species[11,3]=2
species[12,3]=3
#Bismuth and Sulphur
species[13,3]=0
species[14,3]=0
species[15,3]=0
species[16,3]=0
species[17,3]=0
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
species[44,3]=0
species[45,3]=0
species[46,3]=0
#CuBiSO
species[47,3]=2
species[48,3]=7
species[49,3]=1
species[50,3]=1
species[51,3]=2
species[52,3]=1
species[53,3]=1
#BiSO
species[54,3]=0
species[55,3]=0
species[56,3]=0
species[57,3]=0
species[58,3]=0
species[59,3]=0
#CuBiS
species[60,3]=1
species[61,3]=4
species[62,3]=4
#BiCuSO
species[63,3]=1
species[64,3]=3

########### Number of Bismuths Bi #############################################
#Copper
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
species[10,4]=0
species[11,4]=0
species[12,4]=0
#Bismuth
species[13,4]=1
species[14,4]=2
species[15,4]=2
species[16,4]=2
species[17,4]=4
species[18,4]=1
species[19,4]=1
species[20,4]=1
#Sulphur
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
species[44,4]=0
species[45,4]=0
species[46,4]=0
#CuSBiO
species[47,4]=0
species[48,4]=0
species[49,4]=0
species[50,4]=0
species[51,4]=0
species[52,4]=1
species[53,4]=2
#BiSO
species[54,4]=1
species[55,4]=2
species[56,4]=2
species[57,4]=2
species[58,4]=14
species[59,4]=2
#CuBiS
species[60,4]=1
species[61,4]=4
species[62,4]=5
#BiCuSO
species[63,4]=1
species[64,4]=1

########### Number of Sulphurs S  #############################################
#Coppers
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
species[10,5]=0
species[11,5]=0
species[12,5]=0
#Bismuth
species[13,5]=0
species[14,5]=0
species[15,5]=0
species[16,5]=0
species[17,5]=0
species[18,5]=0
species[19,5]=0
species[20,5]=0
#Sulphur
species[21,5]=1
species[22,5]=1
species[23,5]=1
species[24,5]=1
species[25,5]=2
species[26,5]=3
species[27,5]=4
species[28,5]=5
species[29,5]=2
species[30,5]=2
species[31,5]=2
species[32,5]=5
species[33,5]=4
species[34,5]=2
species[35,5]=2
species[36,5]=3
species[37,5]=1
species[38,5]=1
species[39,5]=1
species[40,5]=2
species[41,5]=1
species[42,5]=1
species[43,5]=1
species[44,5]=2
species[45,5]=1
species[46,5]=2
#CuSBiO
species[47,5]=1
species[48,5]=4
species[49,5]=1
species[50,5]=2
species[51,5]=3
species[52,5]=0
species[53,5]=0
#BiSO
species[54,5]=2
species[55,5]=3
species[56,5]=2
species[57,5]=3
species[58,5]=24
species[59,5]=1
#CuBiS
species[60,5]=2
species[61,5]=9
species[62,5]=10
#BiCuSO
species[63,5]=1
species[64,5]=3

######### Number of H2O's #####################################################
#Copper
species[0,6]=0
species[1,6]=1
species[2,6]=1
species[3,6]=0
species[4,6]=0
species[5,6]=4
species[6,6]=2
species[7,6]=3
species[8,6]=1
species[9,6]=2
species[10,6]=1
species[11,6]=2
species[12,6]=4
#Bi
species[13,6]=0
species[14,6]=3
species[15,6]=5
species[16,6]=4
species[17,6]=7
species[18,6]=0
species[19,6]=1
species[20,6]=1
#Sulphur
species[21,6]=0
species[22,6]=0
species[23,6]=0
species[24,6]=0
species[25,6]=0
species[26,6]=0
species[27,6]=0
species[28,6]=0
species[29,6]=3
species[30,6]=3
species[31,6]=3
species[32,6]=6
species[33,6]=6
species[34,6]=4
species[35,6]=4
species[36,6]=6
species[37,6]=3
species[38,6]=3
species[39,6]=3
species[40,6]=6
species[41,6]=4
species[42,6]=4
species[43,6]=4
species[44,6]=8
species[45,6]=5
species[46,6]=5
#CuSBiO
species[47,6]=0
species[48,6]=0
species[49,6]=0
species[50,6]=0
species[51,6]=12
species[52,6]=0
species[53,6]=4
#BiSO
species[54,6]=0
species[55,6]=0
species[56,6]=1
species[57,6]=12
species[58,6]=1
species[59,6]=2
#CuBiS
species[60,6]=0
species[61,6]=0
species[62,6]=0
#BiCuSO
species[63,6]=1
species[64,6]=0

########## Aqueous Ions?????? #################################################
#Copper
species[0,7]=0
species[1,7]=0
species[2,7]=0
species[3,7]=1
species[4,7]=1
species[5,7]=1
species[6,7]=1
species[7,7]=1
species[8,7]=1
species[9,7]=1
species[10,7]=1
species[11,7]=1
species[12,7]=1
#Bismuth
species[13,7]=0
species[14,7]=0
species[15,7]=0
species[16,7]=0
species[17,7]=0
species[18,7]=1
species[19,7]=1
species[20,7]=1
#Sulphur
species[21,7]=0
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
species[44,7]=1
species[45,7]=1
species[46,7]=1
#CuSBiO
species[47,7]=0
species[48,7]=0
species[49,7]=0
species[50,7]=0
species[51,7]=0
species[52,7]=0
species[53,7]=0
#BiSO
species[54,7]=0
species[55,7]=0
species[56,7]=0
species[57,7]=0
species[58,7]=0
species[59,7]=0
#CuBiS
species[60,7]=0
species[61,7]=0
species[62,7]=0
#BiCuSO
species[63,7]=0
species[64,7]=0


#Function to determine species combinations
try:
    combos=load('BiCuOS-speciesCombo.npy')
    num=load('BiCuOS-numberSpecies.npy')
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
    combos=np.zeros((45000,9,3))
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

    save('BiCuOS-speciesCombo.npy', combos)
    save('BiCuOS-numberSpecies.npy', asarray([[combo_num]]))
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
plt.title('Bi-Cu-S Pourbaix Diagram, $\eta_{Bi,Cu,S}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Cu:' +str(composition[1])+'Bi:'+str(composition[2])+'S')


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
plt.title('Bi-Cu-S Pourbaix Diagram, $\eta_{Bi,Cu,S}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Cu:' +str(composition[1])+'Bi:'+str(composition[2])+'S')
chartBox=ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.5, chartBox.height*1.5])
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 0.9), ncol=1)
plt.show()



print('End of Script')