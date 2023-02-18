import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from scipy.stats import linregress

# Constants
g = 9.81

steelElasticMod = 192.73856228694635e9
steelPoisson = 0.35225704496405114
steelYield = 540.981570199884e6
steelShearMod = steelElasticMod/(2+2*steelPoisson)
steelHardConst = 1 # add real numbers
steelHardCoeff = 1

aluminumElasticMod = 73.25376489150618e9
aluminumPoisson = 0.2769261455092313
aluminumYield = 298.6339832889318e6
aluminumShearMod = aluminumElasticMod/(2+2*aluminumPoisson)
aluminumHardConst = 1 # add real numbers
aluminumlHardCoeff = 1

# Dimensional Data
length = 180e-3 # rod length from grip to grip [m]
rodDia = 4.76e-3 # [m]
gripeDia = 52.3e-3 # [m]


# Loading Data
Files = [x for x in listdir('TorsionLab') if '.csv' in x]
Steel = [x for x in listdir('TorsionLab') if '1018' in x]
Aluminum = [x for x in listdir('TorsionLab') if '6061' in x]
Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv('TorsionLab/' + File)

# Minor Calculations
for File in Files:
    Data[File]['Force (N)'] = Data[File]['Force  (kgf)']*g #this adds force in Newtons to the data
    Data[File]['Torque (N.m)'] = Data[File]['Force (N)']*(213/85)*gripeDia #this adds applied torque to the data
    Data[File]['Angle (rad)'] = Data[File]['Angle (deg) ']*(np.pi/180) #this adds applied torque to the data

# Torque Vs Twist Graph
fig = plt.figure(1)
ax = fig.gca()
for File in Aluminum:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4]+' torque') #the label corresponds to what the legend will output
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Aluminum Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

fig = plt.figure(2)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4]+' torque') #the label corresponds to what the legend will output
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Steel Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()


# Shear modulus calculation?
# Ask cherwyn about wether we need to fit our shear modulus or if we can just relate our previously calculated elastic modulus to it


# Calculate yield radius with respect to angle
for File in Files:
    if '1018' in File:
        Data[File]['Yield Radii (m)'] = length*(steelYield/2)/(steelShearMod*Data[File]['Angle (rad)'])
    else:
        Data[File]['Yield Radii (m)'] = length*(aluminumYield/2)/(aluminumShearMod*Data[File]['Angle (rad)'])

# Plot yield radius vs twist angle
fig = plt.figure(3)
ax = fig.gca()
for File in Aluminum:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Yield Radii (m)'],label=File[:len(File)-4]+' Yield radii') #the label corresponds to what the legend will output
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Aluminum Yield Radii vs Twist Angle")
plt.ylabel('Yield Radii (m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

fig = plt.figure(4)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Yield Radii (m)'],label=File[:len(File)-4]+' Yield radii') #the label corresponds to what the legend will output
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Steel Yield Radii vs Twist Angle")
plt.ylabel('Yield Radii (m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

plt.show()
