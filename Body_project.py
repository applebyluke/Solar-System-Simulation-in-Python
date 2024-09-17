# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 20:16:49 2021

@author: apple
"""

import numpy as np

''' A class called Body is used to define any body (star, planet, moon) with given
name, mass, velocity and position.
'''

class Body(object):
    
    def __init__(self, name, mass, velocity, position, initial=False, start=False,
                 load=False):
        
        # The characteristics of the body
        self.name = name
        self.mass = mass
        self.velocity = velocity
        self.position = position
        
        # Set initial conditions for orb_period function
        self.initial = initial
        self.start = start
        self.load = load
           
    def move(self, t, a_d, a):
        # Beeman numerical integration move function by timestamp t 
        fun = self.velocity*t + ((4*a - a_d)*(t**2) / 6)
        self.position = np.add(self.position, fun)
        return self
    
    def v_update(self, t, a_d, a, a_1):
        # Beeman numerical integration method to update the velocity of the body
        fun = (2*a_1 + 5*a - a_d)*t / 6
        self.velocity = np.add(self.velocity, fun)
        return self
    
    def orb_period(self):
        # Function to determine the orbital period of a body (works in conjunction with run_sim)
        ans = 'none'
        
        # If the body is in the first qudrant set initial to true
        if self.position[0] > 0 and self.position[1] > 0:
            self.initial = True
            
        # Once the body crosses the x axis then return 'Start'
        if self.initial == True and self.position[1] < 0:
            self.start = True
            ans = 'Start'
            
        # Once the body crosses back over the x axis set load to be true
        if self.start == True and self.position[1] > 0:
            self.load = True
            
        # Once the body has returned to the start poistion return 'Stop' to stop the timer
        if self.load == True and self.position[1] < 0:
            ans = 'Stop'
        if ans == 'Start' or ans == 'Stop':
            return ans
  
    


