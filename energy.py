import matplotlib
matplotlib.use('TKAgg')
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.stats import sem


def delta_energy_glauber(spin_array, new_spin,itrial,jtrial):
    # spins of four nearest neighbours
    
    #checking boundary conditions
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
      
        
    #calculating energy of  initial spin
    E_final = -J* (new_spin * spin_array[itrial_up,jtrial] + new_spin*spin_array[itrial_down,jtrial] + new_spin*spin_array[itrial,jtrial_left] + new_spin * spin_array[itrial,jtrial_right]) 
    
    old_spin = -1 * new_spin
    
    #calculating energy of final spin
    E_initial = -J*(old_spin * spin_array[itrial_up,jtrial] + old_spin*spin_array[itrial_down,jtrial] + old_spin*spin_array[itrial,jtrial_left] + old_spin * spin_array[itrial,jtrial_right]) 
    E_delta = E_final - E_initial
    return E_delta 


def energy_glauber(spin_array, new_spin,itrial,jtrial):
    #this function computes the energy of a single spin
    
    #checking for initial conditions
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
num_measure = int((nstep-100)/10)




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
energies = []
energies_square = []

energy_error =[]
heat_cap_errors = []

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


#update loop here - for Glauber dynamics
#we are going over all elements in the 2d array n times
if dynamic == "G":
    for kT in kTs:
        energy = []
        energy_square =[]
        for n in range(nstep + 1):
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
                
            # getting rid of initial conditions and autocorrection error between measuerements
            if n >= 100 :
                if(n%10==0):
                    Es = []    #stores the enrgy of a single spin
                    
                    #going over every spin in the lattice to calculate the energy of each spin
                    for ii in range(lx):
                        for jj in range(ly):
                            E = energy_glauber(spin,spin[ii,jj],ii,jj) # more efficient code (np.roll) easy just shift the spin array to the left, right ,down and up and add them to find the change in energy. different function.
                            Es.append(E)
                    
                    #energy of a lattice
                    energy.append((sum(Es)/2))
                    energy_square.append((sum(Es)/2)**2)
                    
        #average energy of a lattice at a temperature.
        energies.append(sum(energy)/len(energy))
        energies_square.append(sum(energy_square)/len(energy_square))
        #error in average energy of lattice at a given tempearture.
        energy_error.append(sem(energy))   
        

        
        # using bootstrap method to get error in heat capacity.
        bootstrap_heat =[]
        for k in range (100):
            
            av_ener_sample = 0 
            av_ener_square_sample = 0
            #selecting 100 random energy values to compute the specific heat capacity
            for n in range(num_measure): # maybe just use 100
                r = np.random.random()
                chosen = int(r* num_measure) 
                av_ener_sample += energy[chosen]
                av_ener_square_sample += energy_square[chosen]
            
            
            av_ener_sample = av_ener_sample/num_measure
            av_ener_square_sample = av_ener_square_sample/num_measure
            #computed suscetibility using the samples magnetisation values
            bootstrap_heat.append((1/((lx*ly)* (kT**2))) * (av_ener_square_sample - av_ener_sample**2) )
        heat_cap_error = (np.mean((np.square(bootstrap_heat))) - (np.mean(bootstrap_heat))**2   )**(1/2)
        heat_cap_errors.append(heat_cap_error)                   
                    
                
    
    
    
    
    
    energies= np.array(energies)
    energies_square = np.array(energies_square)
    C = (1/((lx*ly)*(kTs**2))) *(np.array(energies_square - np.square(energies)))
    critical_capa = np.argmax(C)
    print("Critical temperature is : " + str(kTs[critical_capa]))
    plt.errorbar(kTs,C , xerr =None , yerr = heat_cap_errors  )
    plt.title("Specific heat capacity versus temperature of an Ising model using Glauber dynamics")
    plt.xlabel("Temperature")
    plt.ylabel("Specific heat capacity")
    #plt.errorbar(kTs,energies, xerr = None, yerr = energy_error )  
    #plt.title("Energy versus temperature of an Ising model using Glauber dynamics")
    #plt.xlabel("Temperature")
    #plt.ylabel("Enery")
  
    
    f = open("Energy_glauber.txt", "w")
    for e,kT in zip(energies,kTs):
        f.write(str(e)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("heat_capacity_glauber.txt", "w")
    for c,kT in zip(C,kTs):
        f.write(str(c)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("error_heat_capacity_glauber.txt", "w")
    for heat in (heat_cap_errors ):
        f.write(str(heat) + "\n")
    f.close()
    
    
    
    plt.show()





     
if dynamic == "K":       
    for kT in kTs:
        energy = []
        energy_square = []         
        
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
                    
                    #if the two selected spins are aligned, then there is no change in energy
                    if spin_new_1 == spin_new_2:
                        continue
                    
                    
                    #calculating the change in energy as a result of swapping the spins
                    
                    energy_diff_1 = delta_energy_glauber(spin, spin_new_1,itrial_1,jtrial_1)
                    
                    energy_diff_2 = delta_energy_glauber(spin, spin_new_2,itrial_2,jtrial_2)
                    
                    #accounting for nearest neighbour spins
                    if ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == 1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2)==1) or ((itrial_1 == itrial_2) and abs(jtrial_1 - jtrial_2) == lx -1 ) or ((jtrial_1 == jtrial_2) and abs(itrial_1 - itrial_2) == ly-1):
                        energy_diff_kawa = energy_diff_1 + energy_diff_2 - 4*J*(spin_new_1)*(spin_new_2)
                     
                    #the two selected spins are not neighbours
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

                     
                    
                
            #getting rid of initial conditions and accounting foe auto-correction error in measurements
            if n >= 100 :
                if(n%10==0):
                    Es = []
                    for ii in range(lx):
                        for jj in range(ly):
                            
                            E = energy_glauber(spin,spin[ii,jj],ii,jj)
                            Es.append(E)
                 
                    energy.append((sum(Es)/2))
                    energy_square.append((sum(Es)/2)**2)
                    
    
        energies.append(sum(energy)/len(energy))
        energies_square.append(sum(energy_square)/len(energy_square)) 
        energy_error.append(sem(energy))      



        # using bootstrap method to get error in heat capacity.
        bootstrap_heat =[]
        for k in range (100):
            
            av_ener_sample = 0 
            av_ener_square_sample = 0
            #selecting 100 random energy values to compute the specific heat capacity
            for n in range(num_measure):
                r = np.random.random()
                chosen = int(r* num_measure) 
                av_ener_sample += energy[chosen]
                av_ener_square_sample += energy_square[chosen]
            
            
            av_ener_sample = av_ener_sample/num_measure
            av_ener_square_sample = av_ener_square_sample/num_measure
            #computed suscetibility using the samples magnetisation values
            bootstrap_heat.append((1/((lx*ly)*(kT**2))) * (av_ener_square_sample - av_ener_sample**2) )
        heat_cap_error = (np.mean((np.square(bootstrap_heat))) - (np.mean(bootstrap_heat))**2   )**(1/2)
        heat_cap_errors.append(heat_cap_error)                        
                    
                

    energies= np.array(energies)
    energies_square = np.array(energies_square)
    C = (1/((lx*ly)*(kTs**2))) *(np.array(energies_square - np.square(energies)))
    critical_capa = np.argmax(C)
    print("Critical temperature is : " + str(kTs[critical_capa]))
    plt.errorbar(kTs,C , xerr =None , yerr = heat_cap_errors  )
    plt.title("Specific heat capacity versus temperature of an Ising model using kawasaki dynamics")
    plt.xlabel("Temperature")
    plt.ylabel("Specific heat capacity")
    #plt.errorbar(kTs,energies, xerr = None, yerr = energy_error )  
    #plt.title("Energy versus temperature of an Ising model using kawasaki dynamics")
    #plt.xlabel("Temperature")
    #plt.ylabel("Energy")
  
    
    f = open("Energy_kawasaki.txt", "w")
    for e,kT in zip(energies,kTs):
        f.write(str(e)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("heat_capacity_kawasaki.txt", "w")
    for c,kT in zip(C,kTs):
        f.write(str(c)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("error_heat_capacity_kawasaki.txt", "w")
    for heat in (heat_cap_errors ):
        f.write(str(heat) + "\n")
    f.close()
    
    plt.show()

    
 




