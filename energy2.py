import matplotlib
matplotlib.use('TKAgg')
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation




def delta_energy_glauber(spin_array, new_spin,itrial,jtrial):
    # spins of four nearest neighbours
    #just write a bunh of if statements
    rows,columns = spin_array.shape
   
    
    
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
      
        
    
    E_final = -J* (new_spin * spin_array[itrial_up,jtrial] + new_spin*spin_array[itrial_down,jtrial] + new_spin*spin_array[itrial,jtrial_left] + new_spin * spin_array[itrial,jtrial_right]) 
    
    old_spin = -1 * new_spin
    
    E_initial = -J*(old_spin * spin_array[itrial_up,jtrial] + old_spin*spin_array[itrial_down,jtrial] + old_spin*spin_array[itrial,jtrial_left] + old_spin * spin_array[itrial,jtrial_right]) 
    
    E_delta = E_final - E_initial
    
    return E_delta 


def energy_glauber(spin_array, new_spin,itrial,jtrial):
    # spins of four nearest neighbours
    #just write a bunh of if statements
    rows,columns = spin_array.shape
   
    
    
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
      
        
    
    E_final = -J* (new_spin * spin_array[itrial_up,jtrial] + new_spin*spin_array[itrial_down,jtrial] + new_spin*spin_array[itrial,jtrial_left] + new_spin * spin_array[itrial,jtrial_right]) 
    return E_final


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
energy= []
energy_square = []
energies = []
energies_square = []
kTs = np.linspace(1,3,20)



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
#im=plt.imshow(spin, animated=True)

#update loop here - for Glauber dynamics

#we are going over all elements in the 2d array n times

if dynamic == "K": 
    for kT in kTs:
        energy = []
        energy_square = []         
        for n in range(nstep):
            for i in range(lx):
                for j in range(ly):
                    itrial_1=np.random.randint(0,lx)
                    jtrial_1=np.random.randint(0,ly)
                    spin_new_1=-spin[itrial_1,jtrial_1]
                    
                    itrial_2=np.random.randint(0,lx)
                    jtrial_2=np.random.randint(0,ly)
                    spin_new_2 = -spin[itrial_2,jtrial_2]
                    if spin_new_1 == spin_new_2:
                        continue
                    
                    energy_diff_1 = delta_energy_glauber(spin, spin_new_1,itrial_1,jtrial_1)
                    
                    energy_diff_2 = delta_energy_glauber(spin, spin_new_2,itrial_2,jtrial_2)
                    
                    if ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == 1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2)==1) or ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == lx -1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2) == ly-1):
                        energy_diff_kawa = energy_diff_1 + energy_diff_2 - 4*J*(spin_new_1)*(spin_new_2)
                       
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

                     
                    
                
    
            if n >= 100 :
                if(n%10==0):
                    Es = []
                    for ii in range(lx):
                        for jj in range(ly):
                            
                            E = energy_glauber(spin,spin[ii,jj],ii,jj)
                            
                            #
                            #Es_square.append(E_square)
                            Es.append(E)
        #energy_square.append((sum(Es_square)/2)/len(Es_square))            
                    energy.append((sum(Es)/2))
                    energy_square.append((sum(Es)/2)**2)
                    
    
        energies.append(sum(energy)/len(energy))
        energies_square.append(sum(energy_square)/len(energy_square))                     
                    
                
    
    #each value in the magnetisation list represents the average magnetisation for that SPECIFIC temperature 
    energies= np.array(energies)
    energies_square = np.array(energies_square)
    C = (1/((lx*ly)*kTs)) *(np.array(energies_square - np.square(energies)))
    critical_capa = np.argmax(C)
    print("Critical temperature is : " + str(kTs[critical_capa]))
    plt.scatter(kTs,C)
    #plt.scatter(kTs,energies)   
    plt.show()






"""
if dynamic == "G":
    for kT in kTs:
        energy_square = []
        energy =[]
        
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
                    if energy_diff < 0 :
                        spin[itrial,jtrial] = spin_new
                    else:
                        r = random.random()
                        bolt_weight = np.exp(-(1/kT)*energy_diff)
                        if r <= bolt_weight:
                            spin[itrial,jtrial] = spin_new
                        else: 
                            pass
                
             #occasionally plot or update measurements, eg every 10 sweeps
            if n >= 100 :
                if(n%10==0):
                    for ii in range(lx):
                        for jj in range(ly):
                            
                            E = energy_glauber(spin,spin[ii,jj],ii,jj)
                            
                            E_square = E**2
                            Es_square.append(E_square)
                            Es.append(E)
        energy_square.append((sum(Es_square))/len(Es_square))            
        energy.append(sum(Es)/len(Es))
    #each value in the magnetisation list represents the average magnetisation for that SPECIFIC temperature 
    energy= np.array(energy)
    C = (1/((lx*ly)*np.square(kTs))) *(np.array(energy_square) - np.square(energy))
    critical_capa = np.argmax(C)
    print("Critical temperature is : " + str(kTs[critical_capa]))
    #plt.scatter(kTs,C)
    plt.scatter(kTs,energy)   
    plt.show()


"""


   








