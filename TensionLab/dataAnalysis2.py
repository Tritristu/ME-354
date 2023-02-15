import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from os import listdir
from numpy import mean,std

from scipy.stats import linregress  # This is a linear regression function built into the Scipy library.

# You can call help(linregress) if you'd like to learn more. Or check out https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html

# Material properties
materialConstants = pd.read_csv('materialConstants.csv')
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

Files = [x for x in listdir('experimentalData') if '.csv' in x]
# Dimensioning
sampleDimensions = pd.read_csv(r'sampleDimensions.csv')

# Test Data
Data = {x:{} for x in Files}

for File in Files:
    Data[File] = pd.read_csv('experimentalData/' + File)

i=0
for File in Files:
    emptyList = [[]]*len(Data[File].index) #This is an empty list with the same length as the data file
    #Note, when you multiply Python lists it just copies the list, e.g. [[1]]*3=[[1],[1],[1]]
    Area = sampleDimensions['Thickness(m)'][i]*sampleDimensions['Width(m)'][i]*(10**6)
    #Here's some example new dictionary calculations
    Data[File]['Strain (mm/mm)'] = Data[File]['Axial Strain (mm/mm)']  # this adds another strain column to the data
    Data[File]['Stress (MPa)'] = Data[File]['Load (N)']/Area #this adds stress to the data
    print(sampleDimensions['Thickness(m)'][i])
    print(sampleDimensions['Width(m)'][i])
    print(File)
    print(Area)

    #You'll need to do calculations here
    Data[File]['Instantaneous Area (mm^2)'] = (1-Data[File]['Transverse Strain (mm/mm)'])**2 * Area #Calculate the instantaneous area using the original dimensions and the transverse strain
    Data[File]['True Stress (MPa)'] = Data[File]['Load (N)']/Data[File]['Instantaneous Area (mm^2)'] #Add in the true stress here
    #Data[File]['True Strain (mm/mm)'] = Data[File]['True Stress (MPa)']/ #Add in the true strain here
    i=i+1

#Stress Vs Strain Graph

fig = plt.figure()
ax = fig.gca()
for File in Files:
    ax.plot(Data[File]['Strain (mm/mm)'],Data[File]['Stress (MPa)'],label=File) #the label corresponds to what the legend will output
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Engineering Stress vs Strain")
plt.ylabel('Stress (MPa)')
plt.xlabel('Strain (mm/mm)')
plt.legend() #this turns the legend on, you can manually change entries using legend(['Sample 1', 'Sample 2',...])
plt.show()