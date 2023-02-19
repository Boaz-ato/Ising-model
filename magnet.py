import matplotlib
matplotlib.use('TKAgg')
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.stats import sem


# This function computes the change in energy for a Glauber dynamic
def delta_energy_glauber(spin_array, new_spin,itrial,jtrial):
    
    
    # spins of four nearest neighbours
    #checking for boundary conditions
    
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
      
        
    #computes energy of spin
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
magnetisation= []
magnetisation_error =[]
magnetisation_square = []
sus_errors = []
num_measure = int((nstep-100)/10)
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



#we are going over all elements in the 2d array n times
if dynamic == "G":
    for kT in kTs:
        Ms = []
        Ms_square = []
        for n in range(nstep):
            for i in range(lx):    # not really needed
                for j in range(ly): # not really needed
        
        
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
                
            #getting rid of initial condition (taking measurements after 100 sweeps at or near equilibrium)
            if n >= 100 :
                #getting rid of auto-correction error  (taking measurements at every tenth sweep)
                if(n%10==0):
                    M =abs( spin.sum())            #magnetisation of entire lattice
                    M_square = (spin.sum())**2
                    Ms_square.append(M_square)
                    Ms.append(M)                            # an array of the magnetisations of the lattice
        magnetisation_square.append(sum(Ms_square)/len(Ms_square))            
        magnetisation.append(sum(Ms)/len(Ms))         #an array containing the average of the magnetisation of a lattice at a specific temperature.
        magnetisation_error.append(sem(Ms))         #computes standard error in the values of the magnetisation
        
        
        # using bootstrap method to get error in suscetibility.
        bootstrap_sus =[]
        for k in range (100):
            
            av_mag_sample = 0 
            av_mag_square_sample = 0
            #selecting 100 random magnetisation values to compute the suscetibility 
            for n in range(num_measure): # too many just sample 100
                r = np.random.random()
                chosen = int(r* num_measure) 
                av_mag_sample += Ms[chosen]
                av_mag_square_sample += Ms_square[chosen]
            
            
            av_mag_sample = av_mag_sample/num_measure
            av_mag_square_sample = av_mag_square_sample/num_measure
            #computed suscetibility using the samples magnetisation values
            bootstrap_sus.append((1/((lx*ly)*kT)) * (av_mag_square_sample - av_mag_sample**2) )
        sus_error = (np.mean((np.square(bootstrap_sus))) - (np.mean(bootstrap_sus))**2   )**(1/2)
        sus_errors.append(sus_error)
            
            
            
        
    
#each value in the magnetisation list represents the average magnetisation for that SPECIFIC temperature 
    magnetisation = np.array(magnetisation)
    sus = (1/((lx*ly)*kTs)) *(np.array(magnetisation_square - np.square(magnetisation)))
    critical_sus = np.argmax(sus) #critical value is the maximum value in the suscetibility vs temperature graph.
    print("Critical temperature is : " + str(kTs[critical_sus]))
    plt.errorbar(kTs,sus , xerr =None , yerr = sus_errors  )
    plt.title("Susceptibility versus temperature of an Ising model using Glauber dynamics")
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    #plt.errorbar(kTs,magnetisation, xerr = None, yerr = magnetisation_error )  
    #plt.title("Magnetisation versus temperature of an Ising model using Glauber dynamics")
    #plt.xlabel("Temperature")
    #plt.ylabel("Magnetisation")
  
    
    f = open("magentisation_glauber.txt", "w")
    for mag,kT in zip(magnetisation,kTs):
        f.write(str(mag)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("suscetibility_glauber.txt", "w")
    for s,kT in zip(sus,kTs):
        f.write(str(s)+ " "+ str(kT) + "\n")
    f.close()
    
    
    f = open("error_suscetibility_Glauber.txt", "w")
    for error in (sus_errors  ):
        f.write(str(error) + "\n")
    f.close()
    
    
    f = open("error_magnetisation_Glauber.txt", "w")
    for error in (magnetisation_error):
        f.write(str(error) + "\n")
    f.close()
    
    
    plt.show()





     
if dynamic == "K":        
    for kT in kTs:
        Ms = []
        Ms_square = []     
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
                            
                #getting rid of initial condition 
            if n >= 100 :
                #getting rid of auto-correction error
                if(n%10==0):
                    M =abs( spin.sum())            #magnetisation of entire lattice
                    M_square = (spin.sum())**2
                    Ms_square.append(M_square)
                    Ms.append(M)                            # an array of the magnetisations of the lattice
        magnetisation_square.append(sum(Ms_square)/len(Ms_square))            
        magnetisation.append(sum(Ms)/len(Ms))         #an array containing the average of the magnetisation of a lattice at a specific temperature.
        magnetisation_error.append(sem(Ms))         #computes standard error in the values of the magnetisation
            
            
        # using bootstrap method to get error in suscetibility.
        bootstrap_sus =[]
        for k in range (100):
            
            av_mag_sample = 0 
            av_mag_square_sample = 0
            #selecting 100 random magnetisation values to compute the suscetibility 
            for n in range(num_measure):
                r = np.random.random()
                chosen = int(r* num_measure) 
                av_mag_sample += Ms[chosen]
                av_mag_square_sample += Ms_square[chosen]
            
            
            av_mag_sample = av_mag_sample/num_measure
            av_mag_square_sample = av_mag_square_sample/num_measure
            #computed suscetibility using the samples magnetisation values
            bootstrap_sus.append((1/((lx*ly)*kT)) * (av_mag_square_sample - av_mag_sample**2) )
        sus_error = (np.mean((np.square(bootstrap_sus))) - (np.mean(bootstrap_sus))**2   )**(1/2)
        sus_errors.append(sus_error)
            
            
            
        
    
#each value in the magnetisation list represents the average magnetisation for that SPECIFIC temperature 
    magnetisation = np.array(magnetisation)
    sus = (1/((lx*ly)*kTs)) *(np.array(magnetisation_square - np.square(magnetisation)))
    critical_sus = np.argmax(sus) #critical value is the maximum value in the suscetibility vs temperature graph.
    print("Critical temperature is : " + str(kTs[critical_sus]))
    plt.errorbar(kTs,sus , xerr =None , yerr = sus_errors  )
    plt.title("Susceptibility versus temperature of an Ising model using Kawasaki dynamics")
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    #plt.errorbar(kTs,magnetisation, xerr = None, yerr = magnetisation_error )  
    #plt.title("Magnetisation versus temperature of an Ising model using Kawasaki dynamics")
    #plt.xlabel("Temperature")
    #plt.ylabel("Magnetisation")
  
    
    f = open("magentisation_kawasaki.txt", "w")
    for mag,kT in zip(magnetisation,kTs):
        f.write(str(mag)+ " "+ str(kT) + "\n")
    f.close()
    
    f = open("suscetibility_kawasaki.txt", "w")
    for s,kT in zip(sus,kTs):
        f.write(str(s)+ " "+ str(kT) + "\n")
    f.close()
    
    
    f = open("error_suscetibility_kawasaki.txt", "w")
    for error in (sus_errors  ):
        f.write(str(error) + "\n")
    f.close()
    
    
    f = open("error_magnetisation_kawasaki.txt", "w")
    for error in (magnetisation_error):
        f.write(str(error) + "\n")
    f.close()
    
    plt.show()                  
              



