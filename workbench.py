import numpy as np
from abc import ABC
import math
import pandas as pd
import openpyxl

#Messing around with Arrays
"""
array = np.array([[0,1],[1,2],[2,3],[3,4]])
norms = []
for n in range(len(array)):
    norms.append(np.linalg.norm(array[n]))

print(norms)
"""

#Messing around with Doomsday and planetary angles
"""
ax = (np.array([0,1]))
v = (np.array([2,2]))


normv = v / np.linalg.norm(v)
print(v)
print(normv)
print(ax)
print(np.arccos(np.dot(ax,normv)))
"""

#Messing around with Inheritance
"""
class StemStudent(ABC):
    def __init__(self,name,gpa):
        self.name = name
        self.gpa = gpa
    
    def talk(self):
        pass

    def brag(self):
        print(f'{self.name} says: My gpa is {self.gpa}')

class BioStudent(StemStudent):
    def __init__(self, name, gpa):
        super().__init__(name, gpa)
    
    def talk(self):
        print(f'{self.name} says: I study microbes and shit')
    
    def brag(self):
        return super().brag()

class CompSciStudent(StemStudent):
    def __init__(self, name, gpa):
        super().__init__(name, gpa)
    
    def talk(self):
        print(f'{self.name} says: I study how computers work and shit')
    
    def brag(self):
        return super().brag()
    

class VetStudent(StemStudent):
    def __init__(self, name, gpa,placement):
        super().__init__(name, gpa)
        self.placement = placement
    
    def talk(self):
        print(f'{self.name} says: I study animals and shit')
    
    def brag(self):
        return super().brag()
    
    def where(self):
        print(f"I'll be in {self.placement} in Spring")


class main():

    gabriel = BioStudent('Gabriel','78')
    mandeep = CompSciStudent('Mandeep','probably quite high')
    taiga = VetStudent('Taiga','probably also quite high','Japan')
    gabriel.talk()
    mandeep.talk()
    taiga.talk()
    gabriel.brag()
    mandeep.brag()
    taiga.brag()
    taiga.where()
"""

"""
mylist = pd.DataFrame(np.arange(0,1,0.1),[0,3,4,5,6,1,3,4,3,3])

print(mylist)
"""

vel = [2.6,3.7]
print(vel)

varr = [j*3.154e+7/(1.496e+8) for j in vel]

print(varr)
