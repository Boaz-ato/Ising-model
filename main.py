import matplotlib
matplotlib.use('TKAgg')
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation




# This function computes the change in energy for a Glauber dynamic
def delta_energy_glauber(spin_array, new_spin,itrial,jtrial):
    # spins of four nearest neighbours
    
    rows,columns = spin_array.shape
    
    #checking boundary conditions
    
    itrial_up = itrial -1
    if itrial_up == -1:
        itrial_up = rows - 1 
       
        
    itrial_down = itrial + 1
    if itrial_down == rows:
        itrial_down = 0
      
    
    
    jtrial_left = jtrial -1
    if jtrial_left == -1:
        jtrial_left = columns -1
        
        
    jtrial_right = jtrial + 1
    if jtrial_right == columns:
        jtrial_right = 0
      
        
    #calculating energy of spin
    E_final = -J* (new_spin * spin_array[itrial_up,jtrial] + new_spin*spin_array[itrial_down,jtrial] + new_spin*spin_array[itrial,jtrial_left] + new_spin * spin_array[itrial,jtrial_right]) 
    
    old_spin = -1 * new_spin
    
    E_initial = -J*(old_spin * spin_array[itrial_up,jtrial] + old_spin*spin_array[itrial_down,jtrial] + old_spin*spin_array[itrial,jtrial_left] + old_spin * spin_array[itrial,jtrial_right]) 
    
    E_delta = E_final - E_initial
    
    return E_delta 





J=1.0
kB = 1.0
nstep=10000




size = input("Please enter the length of the two dimensional lattice: ")
try:
    size = int(size)

except ValueError:
    while type(size) != int:
        try:
           size = int(size)
        except ValueError:
            size = input("Try again!" +"\n" + " Please enter the length of the two dimensional lattice(Input must ne an integer. e.g. 50): ")
        

T = input("Please enter the temperature of the system: ")
try:
    T = float(T)

except ValueError:
    while type(T) != float:
        try:
            T = float(T)
        except ValueError:        
            T = input("Try again!" + "\n" + "Please enter the temperature of the system: ")
    

dynamic = input("Please enter the dynamics you want to use for the system. Enter G  for the Glauber method or enter K for the Kawasaki method: ")
dynamic = dynamic.capitalize()
while True:
    if dynamic == "K" or dynamic ==  "G":
        break
    else:
        dynamic = input("Please enter the dynamics you want to use for the system. Enter G  for the Glauber method or enter K for the Kawasaki method: ")
        dynamic = dynamic.capitalize()


lx= size
ly=lx
kT = float(T)

#initialise spins randomly
spin=np.zeros((lx,ly),dtype=float)
for i in range(lx):
    for j in range(ly):
        r=random.random()
        if(r<0.5): 
            spin[i,j]=-1
        if(r>=0.5): 
            spin[i,j]=1

fig = plt.figure()
im=plt.imshow(spin, animated=True)


#update loop here - for Glauber dynamics
#we are going over all elements in the 2d array n times
if dynamic == "G":
    for n in range(nstep):
        for i in range(lx):
            for j in range(ly):
    
    
    #select spin randomly
                itrial=np.random.randint(0,lx)
                jtrial=np.random.randint(0,ly)
                spin_new=-spin[itrial,jtrial]
                
    #compute delta E eg via function (account for periodic BC)
    
                energy_diff = delta_energy_glauber(spin, spin_new,itrial,jtrial)

    #perform metropolis test
                if energy_diff <= 0 :
                    spin[itrial,jtrial] = spin_new
                else:
                    r = random.random()
                    bolt_weight = np.exp(-(1/kT)*energy_diff)
                    if r <= bolt_weight:
                        spin[itrial,jtrial] = spin_new
     
        
        #occasionally plot or update measurements, eg every 10 sweeps   
        if(n%10==0): 
            #update measurements
            #dump output
            f=open('spins.dat','w')
            for i in range(lx):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,spin[i,j]))
            f.close()
            #show animation
            plt.cla()
            im=plt.imshow(spin, animated=True, vmin = -1, vmax = 1)
            plt.draw()
            plt.pause(0.0001)  
       
        
      
if dynamic == "K":         
    #looping over every element in the 2d array n times    
    for n in range(nstep):
        for i in range(lx):
            for j in range(ly):
                #selecting first random spin
                itrial_1=np.random.randint(0,lx)
                jtrial_1=np.random.randint(0,ly)
                spin_new_1=-spin[itrial_1,jtrial_1]
                
                #selecting second random spin
                itrial_2=np.random.randint(0,lx)
                jtrial_2=np.random.randint(0,ly)
                spin_new_2 = -spin[itrial_2,jtrial_2]
                
                #if both spins have the same sign, there is no change in energy and switching them will have no effect
                if spin_new_1 == spin_new_2:
                    continue
                
                #I am considering the change in energy as the same change in energy resulting from doing two consecutive single spin flips in the Glauber dynamics
                
                energy_diff_1 = delta_energy_glauber(spin, spin_new_1,itrial_1,jtrial_1)
                
                energy_diff_2 = delta_energy_glauber(spin, spin_new_2,itrial_2,jtrial_2)
                
                #this checks if the two selected spins are nearest neighbours, as there would be overcounting in the energy.
                if ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == 1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2)==1) or ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == lx -1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2) == ly-1):
                    energy_diff_kawa = energy_diff_1 + energy_diff_2 - 4*J*(spin_new_1)*(spin_new_2)
                
                #if the two particles are not nearest neighbours
                else:
                    energy_diff_kawa = energy_diff_1 + energy_diff_2
                
                
                #perform metropolis test
                spin_1 = spin[itrial_1,jtrial_1]
                spin_2 = spin[itrial_2,jtrial_2]
                if energy_diff_kawa <= 0:
                    spin[itrial_1,jtrial_1] =spin_2
                    spin[itrial_2,jtrial_2] =spin_1
                
                else:
                    r = random.random()
                    bolt_weight = np.exp(-(1/kT)*energy_diff_kawa)
                    if r <= bolt_weight:
                        spin[itrial_1,jtrial_1] = spin_2
                        spin[itrial_2,jtrial_2] = spin_1
                        
                  
              
       
    
    #occasionally plot or update measurements, eg every 10 sweeps
        if(n%10==0): 
    #       update measurements
    #       dump output
            f=open('spins.dat','w')
            for i in range(lx):
                for j in range(ly):
                    f.write('%d %d %lf\n'%(i,j,spin[i,j]))
            f.close()
    #       show animation
            plt.cla()
            im=plt.imshow(spin, animated=True, vmin = -1, vmax = 1)
            plt.draw()
            plt.pause(0.0001)
 
