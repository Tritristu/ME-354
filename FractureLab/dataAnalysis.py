import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # use seaborne instead?
from os import listdir
from scipy.stats import linregress

# Constants
Tiyield = 910e6 # Pa
PCyield = 63.3e6
PMMAyield = 40.2e6

#reading data files FOR SANFORD
Files = [x for x in listdir('.') if '.csv' in x]
Ti0 = [x for x in listdir('.') if 'Ti 0deg' in x]
Ti90 = [x for x in listdir('.') if 'Ti 90deg' in x]
PMMA = [x for x in listdir('.') if 'PMMA' in x]
PC = [x for x in listdir('.') if 'PC' in x]

# # reading data files FOR ERIC
# Files = [x for x in listdir('FractureLab') if '.csv' in x]
# Ti0 = [x for x in listdir('FractureLab') if 'Ti 0deg' in x]
# Ti90 = [x for x in listdir('FractureLab') if 'Ti 90deg' in x]
# PMMA = [x for x in listdir('FractureLab') if 'PMMA' in x]
# PC = [x for x in listdir('FractureLab') if 'PC' in x]

Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv(File) #FOR SANFORD
    # Data[File] = pd.read_csv('FractureLab/' + File)

# Finding Sample Elastic slopes
def slopeFitPlas(displacement, force, a, b):
    slope, intercept, R, P, Err = linregress(displacement[a:b], force[a:b])
    yP = 5 + next(i for i, x in enumerate(displacement) if force[i+5] < 0.95*slope*x + intercept)
    yP1 = yP-5
    return yP1,slope, intercept, R

def slopeFitTi(displacement, force, a, b):
    slope, intercept, R, P, Err = linregress(displacement[a:b], force[a:b])
    yP = 50 + next(i for i, x in enumerate(displacement) if force[i+50] < 0.95*slope*x + intercept)
    yP1 = yP-50
    return yP1,slope, intercept, R

# Calculating slopes, Titanium critical loads
for file in Files:
    if 'Ti 0deg' in file:
        Data[file]['critpoint'],Data[file]['slope'],Data[file]['intercept'],Data[file]['rval'] = slopeFitTi(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],100,450)
        Data[file]['Critical Load (N)'] = Data[file]['Load (N)'][Data[file]['critpoint'][0]]
        print(Data[file]['critpoint'][0], file)
        print(Data[file]['Displacement (mm) '][Data[file]['critpoint'][0]])
    elif 'Ti 90deg' in file:
        Data[file]['critpoint'],Data[file]['slope'],Data[file]['intercept'],Data[file]['rval'] = slopeFitTi(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],100,475)
        Data[file]['Critical Load (N)'] = Data[file]['Load (N)'][Data[file]['critpoint'][0]]
    elif 'PC' in file:
        Data[file]['critpoint'],Data[file]['slope'],Data[file]['intercept'],Data[file]['rval'] = slopeFitPlas(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],100,700)
    else:
        Data[file]['critpoint'],Data[file]['slope'],Data[file]['intercept'],Data[file]['rval'] = slopeFitPlas(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],60,135)


# Calculating max loads/critical load for plastic
for file in Files:
    Data[file]['Max Load (N)'] = max(Data[file]['Load (N)'])
    if file in np.concatenate([PC,PMMA]):
        Data[file]['Critical Load (N)'] = Data[file]['Max Load (N)']

# Calculating geometry factors, (conditional) fracture toughnesses
for file in Files:
    width = Data[file]['Width (mm)'][0]*1e-3
    a = Data[file]['Final crack length (mm)'][0]*1e-3
    base = Data[file]['Base (mm)'][0]*1e-3
    relativeCrackSize = a/width
    geometryFunc = (2+relativeCrackSize)*(0.886 + 4.64*relativeCrackSize - 13.32*relativeCrackSize**2 + 14.72*relativeCrackSize**3 - 5.56*relativeCrackSize**4)/((1-relativeCrackSize)**(3/2))
    Data[file]['Conditional Crack Toughness'] = Data[file]['Critical Load (N)'][0]*geometryFunc/(base*np.sqrt(width))
    if Data[file]['Max Load (N)'][0]/Data[file]['Critical Load (N)'][0] < 1.10:
        if file in np.concatenate([Ti0,Ti90]) and (width - a) >= 2.5*(Data[file]['Conditional Crack Toughness'][0]/Tiyield)**2:
            Data[file]['Fracture Toughness'] = Data[file]['Conditional Crack Toughness']
        elif file in PC and (width - a) >= 2.5*(Data[file]['Conditional Crack Toughness'][0]/PCyield)**2:
            Data[file]['Fracture Toughness'] = Data[file]['Conditional Crack Toughness']
        elif file in PMMA and (width - a) >= 2.5*(Data[file]['Conditional Crack Toughness'][0]/PMMAyield)**2:
            Data[file]['Fracture Toughness'] = Data[file]['Conditional Crack Toughness']
        else:
            Data[file]['Fracture Toughness'] = float("nan")
    else:
        Data[file]['Fracture Toughness'] = float("nan")
    #print(file[0:len(file)-4],'Critical Load (N):',Data[file]['Critical Load (N)'][0],'F:',geometryFunc,'Conditional Crack Toughness:',Data[file]['Conditional Crack Toughness'][0],'Fracture Toughness:',Data[file]['Fracture Toughness'][0])

Ti0Load = []
Ti0Crack = []
Ti90Load = []
Ti90Crack = []
PMMALoad = []
PMMACrack = []
PCLoad = []
PCCrack = []

for File in Files:
    if 'Ti 0deg' in File:
        Ti0Load.append(Data[File]['Critical Load (N)'][Data[File]['critpoint'][0]])
        Ti0Crack.append(Data[File]['Conditional Crack Toughness'][0])
    elif 'Ti 90deg' in File:
        Ti90Load.append(Data[File]['Critical Load (N)'][Data[File]['critpoint'][0]])
        Ti90Crack.append(Data[File]['Conditional Crack Toughness'][0])
    elif 'PMMA' in File:
        PMMALoad.append(Data[File]['Critical Load (N)'][Data[File]['critpoint'][0]])
        PMMACrack.append(Data[File]['Conditional Crack Toughness'][0])
    else:
        PCLoad.append(Data[File]['Critical Load (N)'][Data[File]['critpoint'][0]])
        PCCrack.append(Data[File]['Conditional Crack Toughness'][0])

#
# print('Ti 0deg Critical Load [N]:',np.mean(Ti0Load),'std',np.std(Ti0Load))
# print('Ti 0deg Conditional Fracture Toughness [MPa]:',np.mean(Ti0Crack),'std',np.std(Ti0Crack))
# print('Ti 90deg Critical Load [N]:',np.mean(Ti90Load),'std',np.std(Ti90Load))
# print('Ti 90deg Conditional Fracture Toughness [MPa]:',np.mean(Ti90Crack),'std',np.std(Ti90Crack))
# print('PMMA Critical Load [N]:',np.mean(PMMALoad),'std',np.std(PMMALoad))
# print('PMMA Conditional Fracture Toughness [MPa]:',np.mean(PMMACrack),'std',np.std(PMMACrack))
# print('PC Critical Load [N]:',np.mean(PCLoad),'std',np.std(PCLoad))
# print('PC Conditional Fracture Toughness [MPa]:',np.mean(PCCrack),'std',np.std(PCCrack))


# plotting force vs displacement
fig = plt.figure(1)
ax = fig.gca()
File = 'Ti 0deg 3.csv'
critpoint = Data[File]['critpoint'][0]
ax.scatter(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4],marker='.')
displacements = np.linspace(0,max(Data[File]['Displacement (mm) ']),num=len(Data[File]['Displacement (mm) ']))
ax.plot(displacements,0.95*Data[File]['slope']*displacements + Data[File]['intercept'],label=File[:len(File)-4] + ' 95% Slope',color='r')
ax.scatter(Data[File]['Displacement (mm) '][critpoint],Data[File]['Critical Load (N)'][0],label=File[:len(File)-4] + ' Critical Load',marker='*')

# ax.set_xlim(left = 0,right=0.4)
# ax.set_ylim(bottom = 0,top=20000)
plt.title("Titanium 0\u00b0 Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(2)
ax = fig.gca()
File = 'Ti 90deg 2.csv'
critpoint = Data[File]['critpoint'][0]
ax.scatter(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4],marker='.')
displacements = np.linspace(0,max(Data[File]['Displacement (mm) ']),num=len(Data[File]['Displacement (mm) ']))
ax.plot(displacements,0.95*Data[File]['slope']*displacements + Data[File]['intercept'],label=File[:len(File)-4] + ' 95% Slope',color='r')
ax.scatter(Data[File]['Displacement (mm) '][critpoint],Data[File]['Critical Load (N)'][0],label=File[:len(File)-4] + ' Critical Load',marker='*')
# ax.set_xlim(left = 0,right=0.5)
# ax.set_ylim(bottom = 0,top=20000)
plt.title("Titanium 90\u00b0 Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(3)
ax = fig.gca()
File = 'PMMA 1.csv'
critpoint = Data[File]['critpoint'][0]
ax.scatter(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4],marker='.')
displacements = np.linspace(0,max(Data[File]['Displacement (mm) ']),num=len(Data[File]['Displacement (mm) ']))
ax.plot(displacements,0.95*Data[File]['slope']*displacements + Data[File]['intercept'],label=File[:len(File)-4] + ' 95% Slope',color='r')
ax.scatter(Data[File]['Displacement (mm) '][critpoint],Data[File]['Critical Load (N)'][0],label=File[:len(File)-4] + ' Critical Load',marker='*')
# ax.set_xlim(left = 0,right=0.2)
# ax.set_ylim(bottom = 0,top=300)
plt.title("PMMA Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(4)
ax = fig.gca()
File = 'PC 3.csv'
critpoint = Data[File]['critpoint'][0]
ax.scatter(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4],marker='.')
displacements = np.linspace(0,max(Data[File]['Displacement (mm) ']),num=len(Data[File]['Displacement (mm) ']))
ax.plot(displacements,0.95*Data[File]['slope']*displacements + Data[File]['intercept'],label=File[:len(File)-4] + ' 95% Slope',color='r')
ax.scatter(max(Data[File]['Displacement (mm) ']),Data[File]['Critical Load (N)'][0],label=File[:len(File)-4] + ' Critical Load',marker='*')

# ax.set_xlim(left = 0)
# ax.set_ylim(bottom = 0)
plt.title("Polycarb Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()
plt.show()
