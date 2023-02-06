import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from os import listdir

# Material properties
materialConstants = pd.read_csv(r'data/materialConstants.csv')
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

Files = [x for x in listdir(r'data\experimentalData') if '.csv' in x]
# Dimensioning



# Test Data
Data = {x:{} for x in Files}
print(Data)
# for File in Files:
#     Data[File] = pd.read_csv(r'data\experimentalData/',File)
















