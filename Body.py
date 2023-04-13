import matplotlib
matplotlib.use('TKAgg')
from abc import ABC, abstractmethod

import math
import numpy as np
from numpy.linalg import norm

'''
Body class: Abstract class that possesses a mass, name, and colour for display in the simulation
'''
class Body(ABC):
 
    def __init__(self, name, mass, colour):

        self.name = name
        # mass in kg
        self.m = mass
        # colour - need to strip trailing line return
        self.c = colour.strip()

    @abstractmethod
    def initialise(self, G, p):
        pass

    def updatePos(self, G, dt):
        # keep old position to check for year
        self.r_old = self.r
        
        # update position first: Beeman
        self.r = self.r + self.v*dt + (4*self.a - self.a_old)*dt*dt/6.0
        
    def updateVel(self, G, dt, p):
        # update velocity second: Beeman
        a_new = self.updateAcc(G, p)
        self.v = self.v + (2*a_new + 5*self.a - self.a_old)*dt/6.0
        # now update acc ready for next iteration
        self.a_old = self.a
        self.a = a_new
        
    def updateAcc(self, G, p):
        # update acc (gravitational force law)
        pos = self.r - p.r
        a = -G*p.m*pos/math.pow(norm(pos),3)
        return a
    
    # determine kinetic energy 
    def kineticEnergy(self):
        # ke in J
        ke = (np.dot(self.v, self.v))*self.m/2
        return ke

'''
Planet class: Body subclass that's starting position is defined by an initial orbital radius
'''
class Planet(Body):

    def __init__(self, name, mass, colour, orbit):
        super().__init__(name, mass, colour)
        # orbital radius in m
        self.orbit = orbit
        # set year to zero
        self.year = 0
    
    def initialise(self, G, p):
        # inital position, initial coords = (orbit radius, 0)
        self.r = np.array([self.orbit, 0])
        # inital velocity, tangential to position
        # speed = sqrt(G*marsmass/r)
        if (self.orbit == 0.0):
            #Velocity of sun 
            self.v = np.array([0, 0])
        else:
            #Velocity of orbiting body
            vel = math.sqrt(G*p.m/self.orbit)
            self.v = np.array([0, vel])
        # intial accelatation, using gravitational force law
        if (self.orbit == 0.0):
            self.a = np.array([0, 0])
        else:
            self.a = self.updateAcc(G, p)
        # set acc_old = acc to start Beeman
        self.a_old = self.a

    def getPos(self):
        return (self.r)

    def updatePos(self, G, dt):
        return super().updatePos(G, dt)

    def updateVel(self, G, dt, p):
        return super().updateVel(G, dt, p)

    def updateAcc(self, G, p):
        return super().updateAcc(G, p)

    def newYear(self):
        # update the year when the planet passes the +x axis
        if (self.r_old[1] < 0.0 and self.r[1] >= 0.0):
            self.year +=1
            return True
        else:
            return False

    def kineticEnergy(self):
        return super().kineticEnergy()

'''
Satellite class: Body subclass that's starting velocity is defined by a parameter, and always begins from earth. 
'''

class Satellite(Body): 

    def __init__(self, name, mass, colour, velocity):
        super().__init__(name, mass, colour)

        #Converts from km to AU and seconds to years
        # 1 AU = 1.496e+8 km
        #sc = (3.154e+7/1.496e+8)
        #self.v = np.array([j*sc for j in velocity])
        self.v = np.array(velocity)

    def initialise(self, G, p):
        pass

        #inital position, initial coords = (Earth starting position, 0) + 0.001 earth units at the self.v angle
        self.r = np.array([1.0001, 0])
        # intial accelatation, using gravitational force law
        self.a = self.updateAcc(G, p)
        # set acc_old = acc to start Beeman
        self.a_old = self.a
    
    def updatePos(self, G, dt):
        return super().updatePos(G, dt)

    def updateVel(self, G, dt, p):
        return super().updateVel(G, dt, p)
    
    def updateAcc(self, G, p):
        return super().updateAcc(G, p)
    
    def kineticEnergy(self):
        return super().kineticEnergy()