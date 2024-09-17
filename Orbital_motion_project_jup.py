# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 18:04:18 2021

@author: apple
"""

G = 6.67408 * 10 ** (-11)

import os
path = os.getcwd()
print ("The current working directory is %s" % path)

import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Body_project import Body


# A class for the simulation and display of multi-body systems    

class System(object):
    
    def __init__ (self, bodys):
        
        # An array of the bodies in the system
        self.bodys = bodys
        
        # Create a list to append each of the forces between each body
        self.force = []
        
        # Calculate the force between each body j and the other i bodys 
        for j in range(len(self.bodys)):
            for i in range(len(self.bodys)):
                if j != i:
                    
                    # Calculate the magnatude of force using the formula F = Gm1m2/r^2
                    self.f = ((G * self.bodys[j].mass * self.bodys[i].mass) /
                    (np.linalg.norm(self.bodys[j].position - self.bodys[i].position))**2)
                    
                    # Multiply the magnatude by the direction vector
                    self.f = self.f * ((self.bodys[i].position - 
                    self.bodys[j].position) / np.linalg.norm(self.bodys[j].position - 
                    self.bodys[i].position))
                                                             
                    # Append the force to the force list
                    self.force.append(self.f)
                    
        # Create a list of the net force on each body        
        self.force_combined = []
        n = len(self.bodys)
        
        # For each body sum the forces which act on them
        for k in range(n):
            self.f_c = 0
            for l in range(n-1):
                self.f_c += self.force[l + (n-1)*k]
                
            # Append the force for body k into the list of net forces
            self.force_combined.append(self.f_c)

        '''Create an array where rows are the force on each body and columns
        are the time intervals f_n-1, f_n and f_n+1. To initialize the system
        we simply make f_n-1 = f_n and leave the calculation of f_n+1 for later'''
        
        # Creat an array with size n, 2, 3 with n being the number of bodys 
        # 2 for x and y and components of force and 3 for f_n-1, f_n and f_n+1
        self.force_array = np.zeros([len(self.bodys), 2, 3])
        self.force_array[:,:,0], self.force_array[:,:,1] = self.force_combined, self.force_combined
        
        # Create a list with the mass of each body
        self.mass_list = []
        for i in range(len(self.bodys)):
            self.mass_list.append(self.bodys[i].mass)
        
        # Divide the force array by the mass list to get the acceleration array 
        self.accel_array = np.zeros([len(self.bodys), 2, 3])
        for i in range(len(self.bodys)):
            self.accel_array[i,:,:] = self.force_array[i,:,:] / self.mass_list[i]
        
        # For use in sattelite launch, set launch to be False
        self.launch = False
        
            
    def force_recalc(self):
        # Function to recalcuate the force and acceleration of each body at a new position
        
        '''
        The code in this instance method is almost identical to the code used in
        the __init__ function so to understand this code read the comments above.
        '''
        
        self.force = []
        for j in range(len(self.bodys)):
            for i in range(len(self.bodys)):
                if j != i:
                    self.f = ((G * self.bodys[j].mass * self.bodys[i].mass) /
                    (np.linalg.norm(self.bodys[j].position - self.bodys[i].position))**2)
                    self.f = self.f * ((self.bodys[i].position - 
                    self.bodys[j].position) / np.linalg.norm(self.bodys[j].position - 
                    self.bodys[i].position))
                    self.force.append(self.f)
        self.force_combined = []
        n = len(self.bodys)
        for k in range(n):
            self.f_c = 0
            for l in range(n-1):
                self.f_c += self.force[l + (n-1)*k]
            self.force_combined.append(self.f_c)
        
        # Assign the forces f_n+1
        self.force_array[:,:,2] = self.force_combined
        for i in range(len(self.bodys)):
            self.accel_array[i,:,:] = self.force_array[i,:,:] / self.mass_list[i]
        return self
      
                           
    def time_stamp(self, t):
        # Funciton to estimate the new position and velocity of each body after
        # A given amount of time
        
        # Move all bodys with Beeman numerical integration 
        for i in range(len(self.bodys)):
            self.bodys[i] = self.bodys[i].move(t, self.accel_array[i,:,0], self.accel_array[i,:,1])
            
        # Recalculate force at new position    
        self.force_recalc()
        
        # Calculate velocity of each body 
        for i in range(len(self.bodys)):
            self.bodys[i] = self.bodys[i].v_update(t, self.accel_array[i,:,0],
                                self.accel_array[i,:,1], self.accel_array[i,:,2])
         
        # Set f_n = f_n-1 and f_n+1 = f_n for use in the next time stamp
        self.force_array[:,:,0], self.force_array[:,:,1] = self.force_array[:,:,1], self.force_array[:,:,2]
        
        for i in range(len(self.bodys)):
            self.accel_array[i,:,:] = self.force_array[i,:,:] / self.mass_list[i]
        return self                                              
    
    
    def kinetic(self):
        # Function to calculate the kinetic energy of the system
        
        # For each body i use E = 1/2mv^2 to calculate kinetic energy and add to the total
        E = 0
        for i in range(len(self.bodys)):
            E += self.bodys[i].mass * (np.linalg.norm(self.bodys[i].velocity) ** 2) / 2   
        return E
      
        
    def potential(self):
        # Function to calculate the potential energy of the system
        
        # For each possible pair of bodys in the system the formula for potential energy is used
        U = 0
        for j in range(len(self.bodys)):
            for i in range(len(self.bodys)):
                if j != i:
                    self.pot = ((G * self.bodys[j].mass * self.bodys[i].mass) /
                    (np.linalg.norm(self.bodys[j].position - self.bodys[i].position)))
                    # Divide by two to aviod double counting 
                    U += - self.pot / 2
        return U
    
    
    def energy_total(self):
        # Funciton to calculate the total energy of the system
        
        # Sum the kinetic and potential energy to give the total energy of the system
        return self.kinetic() + self.potential()
         
    
    def sat_launch(self, b0, b1, angle, speed, offset):
        ''' Function to lacunch a sattelite from b0 to b1 when they crate and angle 
        of 'angle' radians with the origin at a given speed. The direction of velocity
        is always the same as b0 unless a negative value is given for speed in which
        case the direction is pi radians from b0's velocity. 
        '''
        
        # Calculate the value of cosine of the angle between b0 and b1
        r = (sum(self.bodys[b0].position * self.bodys[b1].position))
        r /= (np.linalg.norm(self.bodys[b0].position) * np.linalg.norm(self.bodys[b1].position))
        
        det = self.bodys[b0].position[0] * self.bodys[b1].position[1] - (
            self.bodys[b0].position[1] * self.bodys[b1].position[0]) 
        
        # If the angle is within 0.01 of our soecified angle then laucnh the sat
        if (r > np.cos(angle) - 0.001 and r < np.cos(angle) + 0.001 and
        self.launch == False and det < 0):
             
            # Offset the sat by a given amount from b0's to avoid dividing by zero
            offset = offset * (self.bodys[b0].position / np.linalg.norm(self.bodys[b0].position))
           
            # Set the direction of the velocity of the sat 
            direction = self.bodys[b0].velocity / np.linalg.norm(self.bodys[b0].velocity)
            
            # Create an body for the sat
            mars_sat = Body('mars_sat', 100, direction * speed + self.bodys[b0].velocity,
                            self.bodys[b0].position + offset)
            
            # Add the body to the array of bodys 
            self.bodys = np.concatenate((self.bodys, np.array([mars_sat])))
            
            # Call __init__
            self.__init__(self.bodys)
            
            print('Launch!')
            
            # Set self.launch to be True so multiple satellites are not launched
            self.launch = True
        
        return self
            

    def run_sim(self, N, step, orbital=False, energy=False, satellite=False, b0 = 3,
                b1 = 4, angle = 0.767945, speed = 3*10**3, offset = 6.571*10**8, delay=2*10**5, sphere=5.78*10**9):
        ''' Function to run the simulation of the system. Optional arguments are used
        to allow calculation of orbital period, the nergy of the system at each time
        interval and to launch a sattelite from earth to mars. 
        '''
        
        # Initialize the array of positions of bodys used for animation 
        self.bodys_pos = []
        
        # Set launch = False and close = False and sat = 1 to avoid indexing errors in orbital
        sat = 0
        if satellite == True:
            launch = False
            close = False
            sat = 1 
        
        # The file containing energy at each time step is created
        # Likwise the X and Y arrays for the total energy graph are initialized
        if energy == True:
            energyData = open("energy.csv", "w")
            self.X = np.linspace(0, N-1, N)
            
            # Convert X from timesteps to hours 
            self.X *= step / 3600
            self.Y = []  
            
            
        # The inital positions of each body are added to the array used for animation
        for l in range(len(self.bodys)):
            self.bodys_pos.append(np.reshape(self.bodys[l].position, (2,1)))
            
            
        # The arrays for start and end time of the orbital motion calculations are initialized
        if orbital == True:
            T0, T1 = np.zeros(len(self.bodys) + sat), np.zeros(len(self.bodys) + sat)
            
        # Start sat animation at position of Earth
        if satellite == True:
            self.bodys_pos.append(np.reshape(self.bodys[b0].position, (2,1)))
        
        
        # The simulation is ran through N 'step' representing number of seconds 
        for j in range(N):
            self.time_stamp(step)
            
            
            # At each step add the total energy to an array 
            if energy == True:    
                self.Y.append(self.energy_total())
                
            # At each time step the satellite is on standby to launch when the window apears
            if satellite == True and j * step > delay:
                
                self.sat_launch(b0, b1, angle, speed, offset)
                
                # Check if the satellite has launched 
                if self.launch ==  True and launch == False:
                    launch = True
                    
                    # Record the time at which the satellite is launched
                    sat_t0 = j * step 
                                     
                # Record the time at which the satellite gets close to mars
                if (np.linalg.norm(self.bodys[b1].position - self.bodys[-1].position) < sphere
                    and self.launch == True and close == False):
                       
                    close = True
                    sat_t1 = j * step 
                    # Calculate total time of flight in days
                    sat_t = (sat_t1 - sat_t0)
                    sat_t /= 86400
                    
                    print(f'Satellite close to {self.bodys[b1].name} {sat_t} days after launch')
                        
            
            # The position of each body is then added to the animation array
            for i in range(len(self.bodys)):
                xy = np.reshape(self.bodys[i].position, (2,1))
                self.bodys_pos[i] = np.concatenate((self.bodys_pos[i], xy), axis = 1)
                
                
                # The position of each body apart from the first is checked at each time step
                if orbital == True:
                    if self.bodys[i].orb_period() == 'Start' and T0[i] == 0 and i != 0:
                        T0[i] = j
                    if self.bodys[i].orb_period() == 'Stop' and T1[i] == 0 and i != 0:
                        T1[i] = j
                        
                    # The time for one full orbit is the diffenerence between the start and stop times
                    T = T1 - T0
                
                    # Timestamps are converted to earth years
                    T *= step / (3.154*10**7) 
             
            # Add the position of the satellite to the array  
            if satellite == True and launch == False:
                xy = np.reshape(self.bodys[b0].position, (2,1))
                self.bodys_pos[-1] = np.concatenate((self.bodys_pos[-1], xy), axis = 1)                 
        
        # The bodys which have completed at least on full orbit during the simulation have their period printed
        if orbital == True:
            for i in range(len(self.bodys)):
                
                # Bodies that have not completed one full orbit will have T = -n or T = 0
                if T[i] > 0:
                    print(f'The orbital period of {self.bodys[i].name} is {T[i]} Earth years')
        
        # Add the list of total energies at each step to energy.csv
        if energy == True:
            
            # Convert Y to a string
            Y_string = str(self.Y)
            
            # Remove left and right square brackets and add a traling newline
            Y_string = Y_string[1:-1] + '\n'
            
            # Write Y_string to csv
            energyData.write(Y_string)
            
            # Convert X to a list
            X_list = list(self.X)
            
            # Append X as a new row to energy.csv 
            with open('energy.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow(X_list)

    
    def energy_graph(self):
        
        # Take the X and Y from run_sim and plot them on a graph 
        plt.plot(self.X, self.Y)
        plt.xlabel('Hours from start of simulation')
        plt.ylabel('Total energy of system (J)')
        plt.show()
    
    
    def animate(self, i):
        
        # Update position of circles
        for j in range(len(self.bodys)):
            self.patches[j].center = (self.bodys_pos[j][0,i], self.bodys_pos[j][1,i])
            
        return self.patches
                
    
    def display(self, n_frames, sizes, d):
        
        fig = plt.figure()
        ax = plt.axes()
        
        # Create a list of colours so each body can be distinguished 
        colours = ['r', 'g', 'b', 'k', 'c', 'm', 'y']
        
        # Creat a patch for each body 
        self.patches = []
        for k in range(len(self.bodys)):
            self.patches.append(plt.Circle(((self.bodys_pos[k][0,0]), (self.bodys_pos[k][1,0])), sizes[k], color = colours[k], animated = True))
                             
        # Add these patches to the plot
        for i in range(0, len(self.patches)):
            ax.add_patch(self.patches[i])
            
        self.anim = FuncAnimation(fig, self.animate, frames = n_frames, repeat = False, interval = 2, blit = True)
        
        # Set the x and y limits of the plot
        d
        ax.set_xlim(-d, d)
        ax.set_ylim(-d, d)
       
        plt.show()



def sun_main(jupiter=False):
    
    # Opening and reading the csv file that contains the variables for our solar system 
    solarData = open('Solar_System_variables.csv', 'r')
    myReader = csv.reader(solarData)
    
    '''To initialize our solar system we use the mass speed and distance varibales
    from our file to create instances of our Body class and adding them to a list.
    This is done by expressing velocity and position as polar coordinates with 
    randomised theta. This allows us to test many different starting configurations
    of our solar system. The list is then converted to an array.
    '''
    
    solar = []
    for row in myReader:
        theta = 2*np.pi*np.random.rand()
        row[0] = Body((row[0]),
                      float(row[1]),
                      float(row[2]) * np.array([np.cos(theta), np.sin(theta)]),
                      float(row[3]) * np.array([-np.sin(theta), np.cos(theta)]))
        solar.append(row[0])

    if jupiter == False:    
        solar = np.array(solar)
    # Create our solar system by inputting the list of bodys into the Solar system class
    Solar_System = System(solar)
    # Run the simulation with 25000 timesteps of 1 hour with all optional argumnets set to True
    if jupiter == False:
        n_frames = 40000
        d = 3*10**11
        Solar_System.run_sim(40000, 3600, orbital=True, energy=True, satellite=True)
        # Display energy graph
        Solar_System.energy_graph()
    # Run the simulation with jupiter sat launch
    elif jupiter == True:
        n_frames = 100000
        d = 1*10**12
        Solar_System.run_sim(100000, 3600, satellite=True, b0=3,
                             b1=5, angle=1.693, speed=8.8*10**3, delay=6.307*10**7, sphere=4.82*10**10)
    
    # Display the animation with patch sizes given by a list 
    Solar_System.display(n_frames, [3*10**9, 6*10**8, 6*10**8, 6*10**8, 6*10**8, 6*10**8, 6*10**8], d)
    
sun_main(jupiter = True)
