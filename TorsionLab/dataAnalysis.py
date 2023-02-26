import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from scipy.stats import linregress

# Constants
g = 9.81

steelElasticMod = 192.73856228694635e9 # maybe pull from a csv?
steelPoisson = 0.35225704496405114
steelYield = 540.981570199884e6
steelHardExpNom = 0.0875
steelHardCoeffNom = 828e6
# steelHardExpNom = # add extra credit numbers
# steelHardCoeffNom =

aluminumElasticMod = 73.25376489150618e9  # maybe pull from a csv?
aluminumPoisson = 0.2769261455092313
aluminumYield = 298.6339832889318e6
aluminumHardExpNom = 0.0628
aluminumHardCoeffNom = 422e6
# aluminumHardExpEx = # add extra credit numbers
# aluminumHardCoeffEx = 

# Dimensional Data
length = 180e-3 # rod length from grip to grip [m]
rodDia = 4.76e-3 # [m]
gripeDia = 52.3e-3 # [m]
momentOfInertia = np.pi*(rodDia**4)/32

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
    Data[File]['Angle (rad)'] = Data[File]['Angle (deg) ']*(np.pi/180) #this adds angles in radians

# Torque Vs Twist Graph
fig = plt.figure(1)
ax = fig.gca()
for File in Aluminum:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4] + ' Torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Aluminum Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

fig = plt.figure(2)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4] + ' Torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Steel Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

# Shear Modulus Calculation
shearModuli = {}
def shearFit(strain, stress, a, b):
    modulus, intercept, R, P, Err = linregress(strain[a:b], stress[a:b]) 

    return modulus

for File in Files:
    shearStress = Data[File]['Torque (N.m)']*rodDia/(2*momentOfInertia)
    shearStrain = Data[File]['Angle (rad)']*rodDia/(2*length)
    shearModuli[File] = shearFit(shearStrain,shearStress,0,20)

avgSteelShearMod = 0
avgAluminumShearMod = 0
for File in Files:
    if '1018' in File:
        avgSteelShearMod += shearModuli[File]
    else:
        avgAluminumShearMod += shearModuli[File]
avgSteelShearMod = avgSteelShearMod/3
avgAluminumShearMod = avgAluminumShearMod/3

# Calculate yield radius with respect to angle
for File in Files:
    if '1018' in File:
        Data[File]['Yield Radii (m)'] = length*(steelYield/2)/(shearModuli[File]*Data[File]['Angle (rad)'])
    else:
        Data[File]['Yield Radii (m)'] = length*(aluminumYield/2)/(shearModuli[File]*Data[File]['Angle (rad)'])

# Plot yield radius vs twist angle
fig = plt.figure(3)
ax = fig.gca()
for File in Aluminum:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Yield Radii (m)'],label=File[:len(File)-4]+' Yield radii') #the label corresponds to what the legend will output
ax.plot(np.linspace(0,180),(rodDia/2)*np.ones(50),label = 'Specimen Radii')
ax.set_xlim(left=0, right=135)
ax.set_ylim(bottom=0, top=4e-3)
plt.title("Aluminum Yield Radii vs Twist Angle")
plt.ylabel('Yield Radii (m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

fig = plt.figure(4)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Yield Radii (m)'],label=File[:len(File)-4]+' Yield radii') #the label corresponds to what the legend will output
ax.plot(np.linspace(0,180),(rodDia/2)*np.ones(50),label = 'Specimen Radii')
ax.set_xlim(left=0, right=135)
ax.set_ylim(bottom=0, top=4e-3)
plt.title("Steel Yield Radii vs Twist Angle")
plt.ylabel('Yield Radii (m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

# Calculating the theoretical angle-torque plot
def theoreticalPlot(endpoint,length,rodDia,yieldStrength,elasticMod,poissons,hardenCoeff,hardenExp):
    shearMod = elasticMod/(2 + 2*poissons)
    thetaTrans = 2*yieldStrength*length/(shearMod*rodDia*np.sqrt(3)) # rad
    elasticRegion = np.linspace(0,thetaTrans,num=200)
    plasticRegion = np.linspace(thetaTrans,endpoint,num=200)
    elasticTorque = np.pi*shearMod*(rodDia**4)*elasticRegion/((2**5)*length)
    rodYield = yieldStrength*length/(2*shearMod*plasticRegion)
    plasticTorque = (np.pi*yieldStrength*rodYield**3)/(2*np.sqrt(3)) + ((2*np.pi*hardenCoeff)/((hardenExp+3)*np.sqrt(3)))*((plasticRegion/(length*np.sqrt(3)))**hardenExp)*((rodDia/2)**(hardenExp+3) - rodYield**(hardenExp+3))
    angles = np.concatenate([elasticRegion,plasticRegion])*(180/np.pi) # put back into degrees
    torque = np.concatenate([elasticTorque,plasticTorque])
    return angles,torque

anglesSteelTh, torqueSteelTh = theoreticalPlot(8500*(np.pi/180),length,rodDia,steelYield,steelElasticMod,steelPoisson,steelHardCoeffNom,steelHardExpNom)
anglesAluminumTh, torqueAluminumTh = theoreticalPlot(15500*(np.pi/180),length,rodDia,aluminumYield,aluminumElasticMod,aluminumPoisson,aluminumHardCoeffNom,aluminumHardExpNom)

fig = plt.figure(5)
ax = fig.gca()
for File in Aluminum:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4]+' Torque')
ax.plot(anglesAluminumTh,torqueAluminumTh,label='Theoretical torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0, top = 7.5)
plt.title("Aluminum Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

fig = plt.figure(6)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4]+' Torque')
ax.plot(anglesSteelTh,torqueSteelTh,label='Theoretical torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0, top = 15)
plt.title("Steel Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

plt.show()
