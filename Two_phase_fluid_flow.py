# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 14:49:00 2020

@author: Dell
"""


###############################################################################
###############################################################################
###############################################################################
"""#######            multiphase fluid flow equation                        ######"""
###############################################################################
###############################################################################
###############################################################################


import numpy as np
import matplotlib.pyplot as plt 
import math
import pandas as pd

print(   "\t ########################################################################### \t"  )
print(   "\t #############                General solution FOR TWO PHASE SYSTEM             ######## \t"  )
print(   "\t #############    Reservoir Simulation implicit solution to 1D            ############### \t"  )
print(   "\t ######## closed boundary or zero flow rate Neumann boundary condition both side   ############ \t"  )
print(   "\t ##################           injection in first block(flow rate id given)  ############### \t"  )
print(   "\t ##################       production in third block (flow rate id given)    ############### \t"  )
print(   "\t ########################################################################### \t"  )

###############################################################################
"""                   importing data from csv file                          """
print("\n")
heterogeneous_reservoir_data = pd.read_csv('./General_reservoir_description.csv')
print(heterogeneous_reservoir_data.count())
print("\n")

#print(heterogeneous_reservoir_data)



Kx = np.array(heterogeneous_reservoir_data['permeability(md)'])
dx = np.array(heterogeneous_reservoir_data['gridblock_lenght(ft)'])
q_w = np.array(heterogeneous_reservoir_data['Water_flow_rate(ft3/d)'])
q_o = np.array(heterogeneous_reservoir_data['Oil_flow_rate(ft3/d)'])


print("\n absolute permeability in md \n",Kx)
print("\n distance in ft \n",dx)
print("\n water injection flow rate in ft3/d \n",q_w)
print("\n  oil production flow rate in ft3/d \n",q_o)
###############################################################################
###############################################################################

""" Relative permeability data Water Saturation intitial and reduced saturation """
  
Swi = 0.2
print("\n\tinitial water saturation in the reservoir =", Swi)
Swr = 0.2
print("\n\tresidual water saturation in the reservoir =", Swr)
Sw = 0.2
print("\n\t Water saturation in the reservoir =", Sw)


###############################################################################
###############################################################################

number_nodes = len(Kx)-2
print("\n\tBlock node in the reservoir is ", str(number_nodes))

P0 = 1000
print("\n\tThe intial pressure of the reservoir is "+ str(P0)+ "psia")

P_left = 0
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_left))

P_right = 0
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_right))


porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))

viscosity_Water = 1
print("\n\tthe viscosity of water value is "+ str(viscosity_Water) + "cp")


viscosity_oil = 1
print("\n\tthe viscosity of water value is " + str(viscosity_oil) + "cp")

area = 10000
print("\n\tCross sectional area of the reservoir " + str(area) + "ft2")
            

compressibility_w = 1*10**(-5)
print("\n\tcompressibility of water in the reservoir is ", str(compressibility_w))

compressibility_o = 1*10**(-5)
print("\n\tcompressibility of oil in the reservoir is ", str(compressibility_o))

Bw = 1
print("\n\t water formation volume factor is " +str(Bw)+ " rb/stb" )
      

Bo = 1
print("\n\t water formation volume factor is " +str(Bo)+ " rb/stb" )


###############################################################################
###############################################################################
#################### final time for simulation is  ############################

t_final = 3
print("\n\t the reservoir simulation time  should be less than in days is  " +  str(t_final) + "days")

#################### time step  ###############################################

dt_increment = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt_increment)+ "day")
inverse_dt =  ( 1 / dt_increment )

###############################################################################
############            pressure and  boundary condition          #############

pressure_previous = np.ones([number_nodes,1])*P0
print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is\n", str(pressure_previous))

###############################################################################

Saturation_previous = np.ones([number_nodes,1])*Sw
print("\n############## saturation distribution ################\n")
print("\nSaturation distribution at day 0 is\n", str(Saturation_previous))


###############################################################################

print("#########################################################")

def permeability(i,Kx,dx):
    perm = ( dx[i-1] + dx[i] )/( dx[i-1]/Kx[i-1] + dx[i]/Kx[i] )
    return perm

def relative_permeability_water(Sw):
    S = ((Sw - Swi)/(1 - Swi - Swr))
    relperm_w = 0.2*((S)**3)
    #print("Relative permeability of water is = ", relperm_w )
    return relperm_w

"""####    relative permeabiity function from Brooks - Corey    ####"""
def relative_permeability_oil(Sw):
    S = ((Sw - Swi)/(1 - Swi-Swr))
    relperm_o = ((1-S)**3)
    #print("Relative permeability of Oil is = ", relperm_o )
    return relperm_o


def Transmissibility_water(Saturation_previous):
    
    transmissibility_w = np.ones([number_nodes+1,1])
    for i in range(2,number_nodes+2): 
        Sw = Saturation_previous[i-2]
        #print(Sw)
        transmissibility_w[i-1] =  (( permeability(i,Kx,dx)*area*relative_permeability_water(Sw) )/( Bw*viscosity_Water*((dx[i-1] + dx[i])/2) ))
#    Kr_Sw[k] = relative_permeability_water(Saturation_previous[0])
#    Kr2_Sw[k] = relative_permeability_water(Saturation_previous[1])
#    Kr3_Sw[k] = relative_permeability_water(Saturation_previous[2])
    if P_left == 0:
        transmissibility_w[0] = 0*transmissibility_w[0]
    else:
        transmissibility_w[0] = transmissibility_w[0]
        
    if P_right == 0:
        transmissibility_w[number_nodes] = 0*transmissibility_w[number_nodes]
    else:
        transmissibility_w[number_nodes] = transmissibility_w[number_nodes]
    #print(" transmissibility_w is ",transmissibility_w)
    return transmissibility_w

    

def Water_Transmissibility_matrix_water(Saturation_previous):
    transmissibility_w = Transmissibility_water(Saturation_previous)
    transmisibility_matrix_w = np.zeros([number_nodes , number_nodes])
    
    for i in range(1,number_nodes,1):
        transmisibility_matrix_w[i][i] = transmissibility_w[i] + transmissibility_w[i+1]  
        transmisibility_matrix_w[i][i-1] = - transmissibility_w[i]
        transmisibility_matrix_w[i-1][i] = - transmissibility_w[i]    
    transmisibility_matrix_w[0][0] =  transmissibility_w[0] +  transmissibility_w[1]
#    print("\n\ntransmisibility_matrix for Water is \n",transmisibility_matrix_w)
    return transmisibility_matrix_w

######################################################################################


def Transmissibility_oil(Saturation_previous):
    
    transmissibility_o = np.ones([number_nodes+1,1])
    for i in range(2,number_nodes+2): 
        Sw = Saturation_previous[i-2]
        transmissibility_o[i-1] =  ( permeability(i,Kx,dx)*area*relative_permeability_oil(Sw) )/( Bo*viscosity_oil*((dx[i-1] + dx[i])/2) )
    #print("Transmissibility",transmissibility_o)
#    Kr_So[k] = relative_permeability_oil(Saturation_previous[0])
#    Kr2_So[k] = relative_permeability_oil(Saturation_previous[1])
#    Kr3_So[k] = relative_permeability_oil(Saturation_previous[2])
    
    if P_left == 0:
        transmissibility_o[0] = 0*transmissibility_o[0]
    else:
        transmissibility_o[0] = transmissibility_o[0]
        
    if P_right == 0:
        transmissibility_o[number_nodes] = 0*transmissibility_o[number_nodes]
    else:
        transmissibility_o[number_nodes] = transmissibility_o[number_nodes]
    #print(" transmissibility_o is ",transmissibility_o)
    return transmissibility_o    


def Oil_Transmissibility_matrix_oil(Saturation_previous):
    transmissibility_o = Transmissibility_oil(Saturation_previous)
    transmisibility_matrix_o = np.zeros([number_nodes , number_nodes])
    
    for i in range(1,number_nodes,1):
        transmisibility_matrix_o[i][i] = transmissibility_o[i] + transmissibility_o[i+1]  
        transmisibility_matrix_o[i][i-1] = - transmissibility_o[i]
        transmisibility_matrix_o[i-1][i] = - transmissibility_o[i]    
    transmisibility_matrix_o[0][0] =  transmissibility_o[0] +  transmissibility_o[1]
#    print("\n\ntransmisibility_matrix for Oil is \n",transmisibility_matrix_o)
    return transmisibility_matrix_o



def Total_transmissibility_matrix(Saturation_previous): 
    transmisibility_matrix_o = Oil_Transmissibility_matrix_oil(Saturation_previous)
    transmisibility_matrix_w = Water_Transmissibility_matrix_water(Saturation_previous)    
    
    Total_transmissibility_matrix = transmisibility_matrix_w + transmisibility_matrix_o
    
    return Total_transmissibility_matrix

###################################################################################
###################################################################################
###################################################################################
#print("\n############## B_matrix or accumulation matrix  ################\n")
###################################################################################
###################################################################################
Total_compressibility = (Sw)*compressibility_w + (1-Sw)*compressibility_o 
#print("\n Total compressibility = ",Total_compressibility)

B_matrix = np.zeros([number_nodes , number_nodes])
B_inverse_matrix = np.zeros([number_nodes , number_nodes])
B_actual_matrix = np.zeros([number_nodes , number_nodes])


B = np.ones([number_nodes ,1])
for i in range (0,number_nodes):
    B[i] = (area*dx[i+1]*porosity)/Bw
    B_matrix[i][i] = B[i]*Total_compressibility
    B_actual_matrix[i][i]  = B[i]/dt_increment
#print("\n B_matrix is\n", str(B_matrix))

B_actual_matrix = np.linalg.inv(B_actual_matrix)

#print("\nB_actual_matrix\n",B_actual_matrix)

    



###################################################################################
#print("\n############## Q_matrix ################\n")

######################################### injection water matrx   
Q_matrix_w = np.zeros([number_nodes , 1])
for i in range (0,number_nodes):
    Q_matrix_w[i] = q_w[i+1]    
print("\n Q_matrix_w  injection flow rate matrix for Water is\n", str(Q_matrix_w))

######################################### injection water matrx
Q_matrix_o = np.zeros([number_nodes , 1])
for i in range (0,number_nodes):
    Q_matrix_o[i] = q_o[i+1]    
#print("\n Q_matrix_o  production flow rate matrix for oil is\n", str(Q_matrix_o))

Total_flow_rate_Matrix = (Bo/Bw)*Q_matrix_o + Q_matrix_w
print("\n Total flow rate matrix \n", Total_flow_rate_Matrix)


#print("\n##################################################################################")
#      
#B1_Sw = np.ones([t_final , 1])
#B2_Sw = np.ones([t_final , 1])
#B3_Sw = np.ones([t_final , 1])
#
#row = np.ones([t_final , 1])
#
#Kr_Sw = np.ones([t_final , 1])
#Kr_So = np.ones([t_final , 1])
#
#Kr2_Sw = np.ones([t_final , 1])
#Kr2_So = np.ones([t_final , 1])
#
#Kr3_Sw = np.ones([t_final , 1])
#Kr3_So = np.ones([t_final , 1])


x = [166.5,499.5,832.5]
print(x)
x = np.array(x)

for k in range(0, t_final, 1):
    plt.figure(1)
    plt.title('Pressure profile over grid block, Running time ' + str(t_final) + 'days' , fontsize=14  )
    plt.xlabel('Grid block distance  ' , fontsize=16 )
    plt.ylabel('pressure profile ' ,  fontsize=16)
    plt.plot( x , pressure_previous)
    plt.show
#    row[k] = pressure_previous[0]
#    B1_Sw[k] = Saturation_previous[0]
#    B2_Sw[k] = Saturation_previous[1]
#    B3_Sw[k] = Saturation_previous[2]
#    #print("\nB_actual_matrix\n",B_actual_matrix)

#    print("#################( RUNNING TIME )  ####### =  " + str(k) + " DAY  #############")
    pressure_previous = np.dot(np.linalg.inv(( 6.33*10**(-3)*Total_transmissibility_matrix(Saturation_previous) + inverse_dt*B_matrix )), (np.dot((inverse_dt*B_matrix ) , pressure_previous) + Total_flow_rate_Matrix))
    minus_water = - Water_Transmissibility_matrix_water(Saturation_previous)

#    print("\n\nminus_water", minus_water)
#    print(minus_water.shape)
#    print(pressure_previous.shape)
    S_Matrix =  np.dot( minus_water , pressure_previous)
#    print(S_Matrix.shape)
    #print("\n S_Matrix is\n = ",S_Matrix)
    S_Matrix_addition = S_Matrix  + Q_matrix_w
#    print(S_Matrix_addition.shape)
    #print("\n S_Matrix_addition is = ", S_Matrix_addition )
    MM = np.dot(B_actual_matrix, S_Matrix_addition)
#    print(MM.shape)
    #print(" MM  = \n", MM )
    Saturation_previous =  Saturation_previous + MM
    Saturation_previous = np.around(Saturation_previous, decimals = 40) 
#    print(" pressure and saturation after   " + str(k) + " iteration  is ") 
#    print(" previous pressure \n",pressure_previous)
#    print("Saturation previous \n", Saturation_previous)
    #print(Saturation_previous.shape)
    
#
plt.savefig('Pressure vs distanceday2000.png', dpi=1200, bbox_inches='tight')

#   
##    


"""
submission = pd.DataFrame({'id': data_new_test_file['id'], 'left ': prediction})
##Visualize the first 5 rows
submission.head()
filename = 'all_class742308.csv'
submission.to_csv(filename,index=False)
print('Saved file: ' + filename)
#
"""   



"""
#print(Kr_Sw)
Kr_Sw = np.array(Kr_Sw)
Kr_So = np.array(Kr_So)

Kr2_Sw = np.array(Kr2_Sw)
Kr2_So = np.array(Kr2_So)

Kr3_Sw = np.array(Kr3_Sw)
Kr3_So = np.array(Kr3_So)

B1_Sw = np.array(B1_Sw)
B2_Sw = np.array(B2_Sw)
B3_Sw = np.array(B3_Sw)


plt.figure(figsize=(12,8))
plt.legend(["liquid density", "gas density"], prop={"size":20})
plt.plot(B3_Sw, Kr_Sw, color = 'r', label= 'Relative perm of water block 1 ')
plt.plot(B3_Sw, Kr_So, color = 'b', label= 'Relative perm of oil block 1 ')
plt.plot(B3_Sw, Kr2_Sw, color = 'c', label= 'Relative perm of water block 2 ')
plt.plot(B3_Sw, Kr2_So, color = 'g', label= 'Relative perm of oil block 2 ')
plt.plot(B3_Sw, Kr3_Sw, color = 'm', label= 'Relative perm of water block 3 ')
plt.plot(B3_Sw, Kr3_So, color = 'k', label= 'Relative perm of oil block 3 ')

plt.title('Saturation and relative permeability graph running time ' + str(t_final) + 'days' , fontsize=23  )
plt.xlabel('Water Saturation  ' , fontsize=16 )
plt.ylabel('relative permeability of oil and water ' ,  fontsize=16)
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15)
plt.legend()
# save the figure
#plt.savefig('saturation-relative-permeability27.png', dpi=1200, bbox_inches='tight')

plt.show()

"""

###############################################################################
##########              DON't DIsturb          ######################
"""
row = np.array(row)


B1_Sw = np.array(B1_Sw)
B2_Sw = np.array(B2_Sw)
B3_Sw = np.array(B3_Sw)

B1_So = np.ones([t_final , 1])
B2_So = np.ones([t_final , 1])
B3_So = np.ones([t_final , 1])

for h in range(0, t_final, 1):
    B1_So[h] = 1 - B1_Sw[h]
    B2_So[h] = 1 - B2_Sw[h]
    B3_So[h] = 1 - B3_Sw[h]



B1_So = np.array(B1_So)
B2_So = np.array(B2_So)
B3_So = np.array(B3_So)

#print(row)
#print(B1_Sw)

    
plt.figure(figsize=(8,6))
plt.legend(["liquid density", "gas density"], prop={"size":20})
plt.plot(row, B1_Sw, color = 'r', label= 'saturation profile block 1 ')
plt.plot(row, B2_Sw, color = 'g', label= 'saturation profile block 2  ')
plt.plot(row, B3_Sw, color = 'b', label= 'saturation profile block 3  ')
plt.plot(row, B1_So, color = 'c', label= ' oil saturation profile block 1  ')
plt.plot(row, B2_So, color = 'y', label= ' oil saturation profile block 2  ')
plt.plot(row, B3_So, color = 'm', label= ' oil saturation profile block 3  ')

plt.title('pressure and saturation graph', fontsize=23  )
plt.xlabel('pressure psi ' , fontsize=23 )
plt.ylabel('Saturation ',  fontsize=23)
plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=15)
plt.legend()
plt.show()
"""




