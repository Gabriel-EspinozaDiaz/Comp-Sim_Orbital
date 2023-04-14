from Body import Planet
from Body import Satellite

import pandas as pd
import openpyxl
import matplotlib
matplotlib.use('TKAgg')

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy.linalg import norm

'''
Base class to run the orbital simulation
'''
class Solar(object):

    def __init__(self,satellite):

        inputdata = []

        filein = open('parameters-solar-nec.txt', "r")
        for line in filein.readlines():
            if (not line.startswith("#")):
                inputdata.append(line)
        filein.close()

        # simulation parameters
        self.niter = int(inputdata[0]) #Iterations
        self.dt = float(inputdata[1]) #Timesteps
        self.G = float(inputdata[2]) #Gravitational constant

        # list for mars and moons
        self.bodies = []

        # rest of input data is mars and moon data in four line 'chunks'
        # first entry must be the last planet
        for i in range(3, len(inputdata)-4, 4):
            name = inputdata[i]
            mass = float(inputdata[i+1])
            orbit = float(inputdata[i+2])
            colour = inputdata[i+3]
            self.bodies.append(Planet(name, mass, colour, orbit))

        self.peak = False
        if satellite:
            self.bodies.append(Satellite('Perseverance',1e-21,'silver',[4.2,4.3]))
            #Determines the starting distance between Mars and the Satellite
            self.marsDist = np.linalg.norm(np.array([1.524,0])-[1.0001, 0])
            
        # set initial positions and velocities relative to sun
        # sum must be first element in bodies list!
        for i in range(0, len(self.bodies)):
            self.bodies[i].initialise(self.G, self.bodies[0])
        
        # cooldown for checking doomsday
        self.dcool = False

        # storage for system energy
        self.eData = []
        self.kData = []
        self.pData = []

        # Constant for converting from earth masses to kg, seconds to years and AU to meters. 
        # Used for converting arbitrary measure (AU^2 yr^-2) to joules (kg m^2 s^-2) in energy readings
        # 1 earth mass = 5.97219e24 kg
        # 1 AU = 1.496e+11 m
        self.c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)

        # get orbital radius of outermost planet to set size of orbiting bodies and of plot
        # hacky - should really check to see which moon is outermost
        for n in range(len(self.bodies)):
            if type(self.bodies[n]) == Planet:
                self.maxP = self.bodies[n]
    
    def init(self):
        """
        initialiser for animator
        """
        return self.patches 

    def animate(self, i):
        """
        Runs all calculations involved in moving planets, while checking for conditions:
        - New year
        - New Earth year
        - Doomsday alignment
        - Mars distance
        """
        # keep track of time in earth years
        time = (i+1)*self.dt
        self.energy(False)

        # update positions
        for j in range(0, len(self.bodies)):
            self.bodies[j].updatePos(self.G, self.dt)
            self.patches[j].center = self.bodies[j].r
            
        # update velocities
        for j in range(0, len(self.bodies)):
            for k in range(0, len(self.bodies)):
                if (j != k):
                    self.bodies[j].updateVel(self.G, self.dt, self.bodies[k])
        
        # check year and print year if new year for any planet
        for j in range(0, len(self.bodies)):
            if type(self.bodies[j]) == Planet:
                if (self.bodies[j].newYear()):
                    print(self.bodies[j].name.strip() + " " + str(self.bodies[j].year) + " years = " + str(time) + " earth years")
               # in new year is earth year, also print total energy and export total energy
                    if (self.bodies[j].name.strip() == 'earth'):
    
                        energy = self.energy(True)*self.c
                        print('Time = ' + str(time) + ' earth years. Total energy = ' + '{:.3e}'.format(energy) + ' J')

        if not self.peak and type(self.bodies[-1])==Satellite:
            n_marDist = np.linalg.norm(self.bodies[4].getPos()-self.bodies[6].getPos())
            if n_marDist < self.marsDist:
                self.marsDist = n_marDist
            else:
                print(f'lowest distance was {self.marsDist} at {time} years')
                self.peak = True

        # checks doomsday alignment and prints year if true
        if self.dcool == False:
            if self.checkDoomsday():
                print(f'Doomsday at {time} earth years')
        else:
            if time % 0.241 == 0:
                #Resets the cooldown timer approx. every Mercury year
                self.dcool = False
        
        # exports energy readings at the specified point 
        if time == 20.0:
            self.export()

        return self.patches

    def energy(self,returning):
        """
        Returns the system's total energy if True, or stores the current potential, kinetic and total energy in the respective lists if False
        """
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j): #Ensures that the gravitational force being calculated is not 
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double counting
        pe = pe / 2
        totEnergy = ke + pe

        if returning:
            return totEnergy
        else:
            c = self.c
            self.kData.append(ke*c)
            self.pData.append(pe*c)
            self.eData.append(totEnergy*c)

    def calcTotalEnergy(self, i):
        """
        calcluates and prints the total energy of the system at the given iteration
        """
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if (k != j):
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        print('Time = ' + str(i) + ' iterations. Total energy = ' + '{:.3e}'.format(totEnergy)) 

    def checkDoomsday(self):
        """
        Prints out  and sets off a cooldown timer when all planets are aligned
        """
        angles = []
        doomsday = True
        for j in range(1, len(self.bodies)):
            if type(self.bodies[j]) == Planet:
                angles.append(self.bodies[j].getPos() / np.linalg.norm(self.bodies[j].getPos()))
        avg = np.mean(angles,axis=0)
        for j in range(len(angles)):
            #Condition based on 5º in radians, 0.0872665
            if np.arccos(np.clip(np.dot(avg,angles[j]),-1.0,1.0)) > 0.0436332:
                doomsday = False
                break
        if doomsday == True:
            self.dcool = True
            return True
        else:
            return False

    def run(self):
        """
        Runs all steps needed for animation by matplotlib 
        """

        # set up the plot components        
        fig = plt.figure(figsize=[10,10])
        ax = plt.axes()

        # create an array for patches (planet and moons)
        self.patches = []

        # get orbital radius of outermost planet to set size of orbiting bodies and of plot
        # hacky - should really check to see which moon is outermost

        maxOrb = math.sqrt(np.dot(self.maxP.r, self.maxP.r))

        # add the planet, moons and satellites to the Axes and patches
        for i in range(0, len(self.bodies)):
            if type(self.bodies[i]) == Planet:
                if (i == 0):
                    self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.05*maxOrb, color = self.bodies[i].c, animated = True)))
                elif (i > 0 and i < 5):
                    self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.02*maxOrb, color = self.bodies[i].c, animated = True)))
                else:
                    self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.04*maxOrb, color = self.bodies[i].c, animated = True)))
            else:
                    self.patches.append(ax.add_patch(plt.Circle(self.bodies[i].r, 0.01*maxOrb, color = self.bodies[i].c, animated = True)))
        
        # set up the axes
        # scale axes so circle looks like a circle and set limits with border b for prettier plot
        b = 1.2
        lim = maxOrb*b
        print(lim)
        ax.axis('scaled')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
                
        anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.niter, repeat = False, interval = 1, blit= True)

        plt.show()

    def export(self):
        """
        Exports total, kinetic and potential energy to energyData.xlsx
        """
        energyReadings = {'Total Energy':self.eData,'Kinetic Energy':self.kData,'Potential Energy':self.pData}
        export = pd.DataFrame(energyReadings, index = np.arange(0,len(self.eData)*self.dt,self.dt))
        with pd.ExcelWriter('energyData.xlsx') as writer:
            export.to_excel(writer, sheet_name='sheet1')

'''
Solar Subclass that ignores gravitational influence of all bodies but the sun
'''
class SolarNullPlanets(Solar):

    def __init__(self, satellite):
        super().__init__(satellite)
    
    def init(self):
        return super().init()
    
    def animate(self, i):
        """
        Runs all calculations involved in moving planets, while checking for conditions:
        - New year
        - New Earth year
        - Doomsday alignment
        The class differs
        """
        # keep track of time in earth years
        time = (i+1)*self.dt
        self.energy(False)

        # update positions
        for j in range(0, len(self.bodies)):
            self.bodies[j].updatePos(self.G, self.dt)
            self.patches[j].center = self.bodies[j].r
            
        # update velocities based on gravitational pull from the sun alone. 
        for j in range(1, len(self.bodies)):
            self.bodies[j].updateVel(self.G, self.dt, self.bodies[0])
        
        # check year and print year if new year for any planet
        for j in range(0, len(self.bodies)):
            if (self.bodies[j].newYear()):
               print(self.bodies[j].name.strip() + " " + str(self.bodies[j].year) + " years = " + str(time) + " earth years")
               
               # in new year is earth year, also print total energy and export total energy
               if (self.bodies[j].name.strip() == 'earth'):
 
                   energy = self.energy(True)*self.c
                   print('Time = ' + str(time) + ' earth years. Total energy = ' + '{:.3e}'.format(energy) + ' J')

        # checks doomsday alignment and prints year if true
        if self.dcool == False:
            if self.checkDoomsday():
                print(f'Doomsday at {time} earth years')
        else:
            if time % 0.241 == 0:
                #Resets the cooldown timer approx. every Mercury year
                self.dcool = False
        
        # exports energy readings at the specified point 
        if time == 20.0:
            self.export()

        return self.patches

    def energy(self, returning):
        return super().energy(returning)
    
    def calcTotalEnergy(self, i):
        return super().calcTotalEnergy(i)
    
    def checkDoomsday(self):
        return super().checkDoomsday()
    
    def run(self):
        return super().run()
    
    def export(self):
        return super().export()

