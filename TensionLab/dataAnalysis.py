import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from os import listdir

# Material properties
materialConstants = pd.read_csv(r'TensionLab\data\materialConstants.csv')
# 6061-T6 Aluminum
elasticModulusAl = materialConstants.at[0,'Elastic Modulus (Pa)']
yieldStressAl = materialConstants.at[0,'Yield Stress (Pa)']
poissonAl = materialConstants.at[0,'Poisson\'s ratio']
densityAl = materialConstants.at[0,'Density (kg/m^3)']

# 1018 Steel
elasticModulusFe = materialConstants.at[1,'Elastic Modulus (Pa)']
yieldStressFe = materialConstants.at[1,'Yield Stress (Pa)']
poissonFe = materialConstants.at[1,'Poisson\'s ratio']
densityFe = materialConstants.at[1,'Density (kg/m^3)']

# Composite Properties
fiberVol = 0.65 # volume fraction
elasticModulusComp1 = materialConstants.at[3,'Elastic Modulus (Pa)']*(1-fiberVol) + materialConstants.at[2,'Elastic Modulus (Pa)']*fiberVol
elasticModulusComp2 = 1/(fiberVol/materialConstants.at[2,'Elastic Modulus (Pa)'] + (1-fiberVol)/materialConstants.at[3,'Elastic Modulus (Pa)'])
strength1 = fiberVol*materialConstants.at[2,'Yield Stress (Pa)'] + (1-fiberVol)*materialConstants.at[3,'Yield Stress (Pa)']
strength2 = materialConstants.at[3,'Yield Stress (Pa)'] # this is just the strength of the matrix
poissonComp12 = fiberVol*materialConstants.at[2,'Poisson\'s ratio'] + (1-fiberVol)*materialConstants.at[3,'Poisson\'s ratio']
poissonComp21 = (elasticModulusComp2/elasticModulusComp1)*poissonComp12

# Exmperimental Data

Files = [x for x in listdir('TensionLab\data\experimentalData') if '.csv' in x]
# Dimensioning
materialConstants = pd.read_csv(r'TensionLab\data\sampleDimensions.csv')
# Test Data
Data = {x:{} for x in Files}

for File in Files:
    Data[File] = pd.read_csv('TensionLab\data\experimentalData\\' + File,sep=";")


for File in Files:
    emptyList = [[]]*len(Data[File].index) #This is an empty list with the same length as the data file
    #Note, when you multiply Python lists it just copies the list, e.g. [[1]]*3=[[1],[1],[1]]
    
    #Here's some example new dictionary calculations
    Data[File]['Stress (MPa)'] = Data[File]['Load (N)']#/Area #this adds stress to the data
    
    #You'll need to do calculations here
    Data[File]['Instantaneous Area (mm^2)'] = emptyList #Calculate the instantaneous area using the original dimensions and the transverse strain
    Data[File]['True Stress (MPa)'] = emptyList #Add in the true stress here
    Data[File]['True Strain (mm/mm)'] = emptyList #Add in the true strain here,




