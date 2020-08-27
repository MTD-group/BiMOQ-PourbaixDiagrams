#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:17:41 2018

@author: laurenwalters
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
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
k_b=8.617333262145*10**-5

#Array showing the composition of Ag:Bi:Se
composition=np.array([10,1,10])

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
##########################PBEsol with SOC######################################
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
Ee_O2=-10.281123;
Ee_Se=-3.82509;
###############################################################################
########### MULTICOMPONENT SPECIES ############################################
#Vibrational Energies in eV/f.u.
#Calculated with PBEsol+SOC
Ee_Ag5BiO4=-46.0621275;
Ee_Ag2BiO3=-30.539865;
Ee_AgBiO3=-26.74457833;

Ee_AgSe=-6.955571;
Ee_Ag2Se=-10.6994305;
Ee_Ag2SeO3=-28.97367;
Ee_Ag2SeO4=-34.3132895;

Ee_Bi3Se4=-32.7784733;
Ee_Bi4Se3=-33.8044333;
Ee_Bi2Se2=-19.06480267;
Ee_BiSe2=-13.56499625;
Ee_Bi2Se3=-23.401805;
Ee_Bi2SeO2=-29.049585;
Ee_Bi2SeO5=-48.87748;
Ee_Bi2SeO3_3=-82.53894;

Ee_AgBiSeO=-19.7504955;
Ee_AgBiSe2=-16.853356;
###############################################################################
###### Vibrational Energy #####################################################
###############################################################################
#From PBEsol Phonon Calculations
Fvib_O2=-0.272;
F_rot_trans_O2=0.099;
Ftot_O2=Fvib_O2+F_rot_trans_O2; 
F_H = .202;

Fvib_Se=-0.045190058131115

Fvib_AgO=  0.01360462162
Fvib_Ag2O=  -0.07478725194
Fvib_Ag2O3= 0.17837389128
Fvib_Ag= -0.038443004392 

Fvib_Bi=-0.0761976993239
Fvib_Bi2O3=-0.057653546889
Fvib_Bi2O5=0.14677315404
Fvib_Bi2O4=0.12231438709
Fvib_Bi4O7=0.08741679245

Fvib_Ag5BiO4=-0.12861399053
Fvib_Ag2BiO3=-0.098643737028
Fvib_AgBiO3=0.0140751179482

Fvib_AgSe=-0.112913176963
Fvib_Ag2Se=-0.17376246767031
Fvib_Ag2SeO3=-0.0023185734781
Fvib_Ag2SeO4=-0.026133015641

Fvib_Bi3Se4=-0.417899867987
Fvib_Bi4Se3=-0.43313969435
Fvib_Bi2Se2=-0.2343788142756
Fvib_BiSe2=-0.1523250019192
Fvib_Bi2Se3=-0.25487956184
Fvib_Bi2SeO2=-0.05124595680
Fvib_Bi2SeO5=0.12343256096
Fvib_Bi2SeO3_3=0.2336827957

Fvib_AgBiSeO=-0.104731934446
Fvib_AgBiSe2=-0.22428689674
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

dGf_Ag5BiO4=(Ee_Ag5BiO4+Fvib_Ag5BiO4) -5*(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-2.0*((Ee_O2)-Ftot_O2);
dGf_Ag2BiO3=(Ee_Ag2BiO3+Fvib_Ag2BiO3) -2*(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-1.5*((Ee_O2)-Ftot_O2);
dGf_AgBiO3=(Ee_AgBiO3+Fvib_AgBiO3)    -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-1.5*((Ee_O2)-Ftot_O2);

dGf_AgSe=(Ee_AgSe+Fvib_AgSe)          -(Ee_Ag+Fvib_Ag)-(Ee_Se+Fvib_Se);
dGf_Ag2Se=(Ee_Ag2Se+Fvib_Ag2Se)       -2*(Ee_Ag+Fvib_Ag)-(Ee_Se+Fvib_Se);
dGf_Ag2SeO3=(Ee_Ag2SeO3+Fvib_Ag2SeO3) -2*(Ee_Ag+Fvib_Ag)-(Ee_Se+Fvib_Se)-1.5*((Ee_O2)-Ftot_O2);
dGf_Ag2SeO4=(Ee_Ag2SeO4+Fvib_Ag2SeO4) -2*(Ee_Ag+Fvib_Ag)-(Ee_Se+Fvib_Se)-2.0*((Ee_O2)-Ftot_O2);

dGf_Bi3Se4=(Ee_Bi3Se4+Fvib_Bi3Se4)          -4*(Ee_Se+Fvib_Se)-3*(Ee_Bi+Fvib_Bi);
dGf_Bi4Se3=(Ee_Bi4Se3+Fvib_Bi4Se3)          -3*(Ee_Se+Fvib_Se)-4*(Ee_Bi+Fvib_Bi);
dGf_Bi2Se2=(Ee_Bi2Se2+Fvib_Bi2Se2)          -2*(Ee_Se+Fvib_Se)-2*(Ee_Bi+Fvib_Bi);
dGf_BiSe2=(Ee_BiSe2+Fvib_BiSe2)             -2*(Ee_Se+Fvib_Se)-(Ee_Bi+Fvib_Bi);
dGf_Bi2Se3=(Ee_Bi2Se3+Fvib_Bi2Se3)          -3*(Ee_Se+Fvib_Se)-2*(Ee_Bi+Fvib_Bi);
dGf_Bi2SeO2=(Ee_Bi2SeO2+Fvib_Bi2SeO2)       -(Ee_Se+Fvib_Se)-2*(Ee_Bi+Fvib_Bi)-1.0*((Ee_O2)-Ftot_O2);
dGf_Bi2SeO5=(Ee_Bi2SeO5+Fvib_Bi2SeO5)       -(Ee_Se+Fvib_Se)-2*(Ee_Bi+Fvib_Bi)-2.5*((Ee_O2)-Ftot_O2);
dGf_Bi2SeO3_3=(Ee_Bi2SeO3_3+Fvib_Bi2SeO3_3) -3*(Ee_Se+Fvib_Se)-2*(Ee_Bi+Fvib_Bi)-4.5*((Ee_O2)-Ftot_O2);

dGf_AgBiSeO=(Ee_AgBiSeO+Fvib_AgBiSeO) -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-(Ee_Se+Fvib_Se)-0.5*((Ee_O2)-Ftot_O2);
dGf_AgBiSe2=(Ee_AgBiSe2+Fvib_AgBiSe2) -(Ee_Ag+Fvib_Ag)-(Ee_Bi+Fvib_Bi)-2*(Ee_Se+Fvib_Se);

#Set the reference values
dGf_Ag=0.0;
dGf_Bi=0.0;
dGf_Se=0.0;
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

#Elemental Selenium Species
dGf_H2Se= 0.7978997
dGf_HSe_Minus=1.01906
dGf_Se_2Minus= 1.8473
dGf_H2SeO3=-4.4144669
dGf_HSeO3_Minus= -4.26269
dGf_SeO3_2Minus=-3.8737164
dGf_H2SeO4=-4.571445
dGf_HSeO4_Minus=-4.692
dGf_SeO4_2Minus=-4.5757814

#Water
dGf_H2O=-2.458
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

dGf_Se=   dGf_Se*F;   
dGf_H2Se= dGf_H2Se*F;          
dGf_HSe_Minus=dGf_HSe_Minus*F;   
dGf_Se_2Minus= dGf_Se_2Minus*F;  
dGf_H2SeO3= dGf_H2SeO3*F;             
dGf_HSeO3_Minus= dGf_HSeO3_Minus*F;  
dGf_SeO3_2Minus= dGf_SeO3_2Minus*F;   
dGf_H2SeO4= dGf_H2SeO4*F;            
dGf_HSeO4_Minus= dGf_HSeO4_Minus*F;  
dGf_SeO4_2Minus= dGf_SeO4_2Minus*F;  

dGf_Ag5BiO4=dGf_Ag5BiO4*F;
dGf_Ag2BiO3=dGf_Ag2BiO3*F;
dGf_AgBiO3=dGf_AgBiO3*F;

dGf_AgSe=dGf_AgSe*F;
dGf_Ag2Se=dGf_Ag2Se*F;
dGf_Ag2SeO3=dGf_Ag2SeO3*F;
dGf_Ag2SeO4=dGf_Ag2SeO4*F;

dGf_Bi3Se4=dGf_Bi3Se4*F;
dGf_Bi4Se3=dGf_Bi4Se3*F;
dGf_Bi2Se2=dGf_Bi2Se2*F;
dGf_BiSe2=dGf_BiSe2*F;
dGf_Bi2Se3=dGf_Bi2Se3*F;
dGf_Bi2SeO2=dGf_Bi2SeO2*F;
dGf_Bi2SeO5=dGf_Bi2SeO5*F;
dGf_Bi2SeO3_3=dGf_Bi2SeO3_3*F;

dGf_AgBiSeO=dGf_AgBiSeO*F;
dGf_AgBiSe2=dGf_AgBiSe2*F;

dGf_H2O=  dGf_H2O*F;                   
###############################################################################
###############################################################################
###############################################################################

###############################################################################
###############   Populate the species matrix  ################################
###############################################################################
species=np.zeros((45,8))

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

species[18,0]=dGf_Se
species[19,0]=dGf_H2Se
species[20,0]=dGf_HSe_Minus
species[21,0]=dGf_Se_2Minus
species[22,0]=dGf_H2SeO3
species[23,0]=dGf_HSeO3_Minus
species[24,0]=dGf_SeO3_2Minus
species[25,0]=dGf_H2SeO4
species[26,0]=dGf_HSeO4_Minus
species[27,0]=dGf_SeO4_2Minus

species[28,0]=dGf_AgSe
species[29,0]=dGf_Ag2Se
species[30,0]=dGf_Ag2SeO3
species[31,0]=dGf_Ag2SeO4
species[32,0]=dGf_AgBiSeO;
species[33,0]=dGf_AgBiSe2;
species[34,0]=dGf_Ag5BiO4
species[35,0]=dGf_Ag2BiO3
species[36,0]=dGf_AgBiO3

species[37,0]=dGf_Bi3Se4
species[38,0]=dGf_Bi4Se3
species[39,0]=dGf_Bi2Se2
species[40,0]=dGf_BiSe2
species[41,0]=dGf_Bi2Se3
species[42,0]=dGf_Bi2SeO2
species[43,0]=dGf_Bi2SeO5
species[44,0]=dGf_Bi2SeO3_3

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
#Selenium
species[18,1]=0
species[19,1]=-2
species[20,1]=-2
species[21,1]=-2
species[22,1]=4
species[23,1]=4
species[24,1]=4
species[25,1]=6
species[26,1]=6
species[27,1]=6
#AgSeBiO
species[28,1]=0
species[29,1]=0
species[30,1]=6
species[31,1]=8
species[32,1]=2
species[33,1]=0
species[34,1]=8
species[35,1]=6
species[36,1]=6
#BiSeO
species[37,1]=0
species[38,1]=0
species[39,1]=0
species[40,1]=0
species[41,1]=0
species[42,1]=4
species[43,1]=10
species[44,1]=18

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
#Selenium
species[18,2]=0
species[19,2]=-2
species[20,2]=-1
species[21,2]=0
species[22,2]=4
species[23,2]=5
species[24,2]=6
species[25,2]=6
species[26,2]=7
species[27,2]=8
#AgSeBiO
species[28,2]=0
species[29,2]=0
species[30,2]=6
species[31,2]=8
species[32,2]=2
species[33,2]=0
species[34,2]=8
species[35,2]=6
species[36,2]=6
#BiSeO
species[37,2]=0
species[38,2]=0
species[39,2]=0
species[40,2]=0
species[41,2]=0
species[42,2]=4
species[43,2]=10
species[44,2]=18

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
#Selenium
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
#AgSeBiO
species[28,3]=1
species[29,3]=2
species[30,3]=2
species[31,3]=2
species[32,3]=1
species[33,3]=1
species[34,3]=5
species[35,3]=2
species[36,3]=1
#BiSeO
species[37,3]=0
species[38,3]=0 
species[39,3]=0
species[40,3]=0
species[41,3]=0
species[42,3]=0
species[43,3]=0
species[44,3]=0

########### Number of Bismuths Bi #############################################
#Silver
species[0,4]=0;
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

#Selenium
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
#AgSeBiO
species[28,4]=0
species[29,4]=0
species[30,4]=0
species[31,4]=0
species[32,4]=1
species[33,4]=1
species[34,4]=1
species[35,4]=1
species[36,4]=1
#BiSeO
species[37,4]=3
species[38,4]=4
species[39,4]=2
species[40,4]=1
species[41,4]=2
species[42,4]=2
species[43,4]=2
species[44,4]=2

########### Number of Selenium Se  ############################################
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
#Selenium
species[18,5]=1
species[19,5]=1
species[20,5]=1
species[21,5]=1
species[22,5]=1
species[23,5]=1
species[24,5]=1
species[25,5]=1
species[26,5]=1
species[27,5]=1
#AgSeBiO
species[28,5]=1
species[29,5]=1
species[30,5]=1
species[31,5]=1
species[32,5]=1
species[33,5]=2
species[34,5]=0
species[35,5]=0
species[36,5]=0
#BiSeO
species[37,5]=4
species[38,5]=3
species[39,5]=2
species[40,5]=2
species[41,5]=3
species[42,5]=1
species[43,5]=1
species[44,5]=3

######### Number of H2O's #####################################################
#Ag
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
#Selenium
species[18,6]=0
species[19,6]=0
species[20,6]=0
species[21,6]=0
species[22,6]=3
species[23,6]=3
species[24,6]=3
species[25,6]=4
species[26,6]=4
species[27,6]=4
#AgSeBiO
species[28,6]=0
species[29,6]=0
species[30,6]=3
species[31,6]=4
species[32,6]=1
species[33,6]=0
species[34,6]=4
species[35,6]=3
species[36,6]=3
#BiSeO
species[37,6]=0
species[38,6]=0
species[39,6]=0
species[40,6]=0
species[41,6]=0
species[42,6]=2
species[43,6]=5
species[44,6]=9

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
#Selenium
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
#AgSeBiO
species[28,7]=0
species[29,7]=0
species[30,7]=0
species[31,7]=0
species[32,7]=0
species[33,7]=0
species[34,7]=0
species[35,7]=0
species[36,7]=0
#SeBiO
species[37,7]=0
species[38,7]=0
species[39,7]=0
species[40,7]=0
species[41,7]=0
species[42,7]=0
species[43,7]=0
species[44,7]=0

#Function to determine species combinations
try:
    combos=load('BiAgOSe-speciesCombo.npy')
    num=load('BiAgOSe-numberSpecies.npy')
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
    combos=np.zeros((20000,9,3))
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

    save('BiAgOSe-speciesCombo.npy', combos)
    save('BiAgOSe-numberSpecies.npy', asarray([[combo_num]]))
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
            current_mu=((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
             
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
            current_mu=current_mu+((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
            
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
            current_mu=current_mu+((current_eng+current_aquI-current_ele-current_H-current_H2O)/current_NumEle);
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
            elif(flag[k,0]==muValues[i,j,1] and flag[k,1]==muValues[i,j,2]and flag[k,2]==muValues[i,j,0]):
                ax.plot(pH,U,'.', color = [flag[k,3],flag[k,4],flag[k,5]],markersize=4)
                #break loop, the color is found
                k=len(flag)+1
                l=1
            elif(flag[k,0]==muValues[i,j,2] and flag[k,1]==muValues[i,j,0]and flag[k,2]==muValues[i,j,1]):
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
plt.title('Bi-Ag-Se Pourbaix Diagram, $\eta_{Bi,Ag,Se}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Ag:' +str(composition[1])+'Bi:'+str(composition[2])+'Se')


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
plt.title('Bi-Ag-Se Pourbaix Diagram, $\eta_{Bi,Ag,Se}=10^{-'+str(eta)+'}$, '+str(composition[0])+'Ag:' +str(composition[1])+'Bi:'+str(composition[2])+'Se')
chartBox=ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.5, chartBox.height*1.5])
ax.legend(loc='upper center', bbox_to_anchor=(1.3, 0.9), ncol=1)
plt.show()



print('End of Script')