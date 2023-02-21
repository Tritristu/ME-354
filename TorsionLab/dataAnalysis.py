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
momentOfInertia = np.pi*(rodDia**4)/64

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
    ax.plot(Data[File]['Angle (rad)'],Data[File]['Torque (N.m)'],label=File[:len(File)-4] + ' Torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Aluminum Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (rad)')
plt.legend()

fig = plt.figure(2)
ax = fig.gca()
for File in Steel:
    ax.plot(Data[File]['Angle (rad)'],Data[File]['Torque (N.m)'],label=File[:len(File)-4] + ' Torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Steel Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (rad)')
plt.legend()

# Shear Modulus Calculation
shearModuli = {}
def shearFit(angle, torque, a, b):
    nu, C, R, P, Err = linregress(torque[a:b], angle[a:b])  # The data outputs the slope (nu), intercept (C), regression (R) value, P-value and standard error

    # Y = [0.0, torque[round(1.5 * b)]]
    # X = [(y - C) / nu for y in Y]  # these are points that you can plot to visualize the data being fit, inverted from y=nu*x+C, x=(y-C)/nu
    return nu

for File in Files:
    shearStress = Data[File]['Torque (N.m)']*rodDia/(2*momentOfInertia)
    shearStrain = Data[File]['Angle (deg) ']*rodDia/(2*length)
    shearModuli[File] = shearFit(shearStress,shearStrain,0,10)

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
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
# ax.set_xlim(left=0, right=0.01)
ax.set_ylim(bottom=0, top=rodDia)
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
# ax.set_xlim(left=0, right=0.01)
ax.set_ylim(bottom=0, top=rodDia)
plt.title("Steel Yield Radii vs Twist Angle")
plt.ylabel('Yield Radii (m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()

# Calculating the theoretical angle-torque plot
def theoreticalPlot(endpoint,yieldStrength,shearMod,hardenCoeff,hardenExp):
    length = 180e-3
    rodDia = 4.76e-3
    thetaTrans = yieldStrength*length/(shearMod*rodDia)
    elasticRegion = np.linspace(0,thetaTrans)
    plasticRegion = np.linspace(thetaTrans,endpoint,num=100)
    elasticTorque = np.pi*shearMod*rodDia*elasticRegion*(np.pi/180)/(64*length)
    rodYield = (yieldStrength/2)*length/(shearMod*plasticRegion*(np.pi/180))
    plasticTorque = ((2*np.pi*hardenCoeff)/((hardenExp+3)*np.sqrt(3)))*((plasticRegion*(np.pi/180)/(length*np.sqrt(3)))**3)*(rodDia**(hardenExp+3)-rodYield**(hardenExp+3))
    angles = np.concatenate([elasticRegion,plasticRegion]) # having trouble concatenating these two lists
    torque = np.concatenate([elasticTorque,plasticTorque])
    return angles,torque

# anglesSteelTh, torqueSteelTh = theoreticalPlot(1000,steelYield,avgSteelShearMod,steelHardCoeff,steelHardExp)
# anglesAluminumTh, torqueAluminumTh = theoreticalPlot(1600,aluminumYield,avgAluminumShearMod,aluminumHardCoeffNom,aluminumHardExpNom)

fig = plt.figure(5)
ax = fig.gca()
# for File in Aluminum:
#     ax.plot(Data[File]['Angle (deg) '],Data[File]['Torque (N.m)'],label=File[:len(File)-4]+' torque')
# ax.plot(anglesAluminumTh,torqueAluminumTh,label='Theortical torque')
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Aluminum Torque vs Twist Angle")
plt.ylabel('Torque (N.m)')
plt.xlabel('Twist Angle (deg)')
plt.legend()


plt.show()
