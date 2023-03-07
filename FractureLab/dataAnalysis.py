import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # use seaborne instead?
from os import listdir
from scipy.stats import linregress

# Constants
Tiyield = 240e6 # Pa
PCyield = 66e6
PMMAyield = 70e6

# reading data files
Files = [x for x in listdir('FractureLab') if '.csv' in x]
Ti0 = [x for x in listdir('FractureLab') if 'Ti 0deg' in x]
Ti90 = [x for x in listdir('FractureLab') if 'Ti 90deg' in x]
PMMA = [x for x in listdir('FractureLab') if 'PMMA' in x]
PC = [x for x in listdir('FractureLab') if 'PC' in x]

Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv('FractureLab/' + File)

# plotting force vs displacement
fig = plt.figure(1)
ax = fig.gca()
for File in Ti0:
    ax.plot(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4])
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Titanium 0\u00b0 Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(2)
ax = fig.gca()
for File in Ti90:
    ax.plot(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4])
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Titanium 90\u00b0 Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(3)
ax = fig.gca()
for File in PMMA:
    ax.plot(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4])
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("PMMA Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()

fig = plt.figure(4)
ax = fig.gca()
for File in PC:
    ax.plot(Data[File]['Displacement (mm) '],Data[File]['Load (N)'],label=File[:len(File)-4])
ax.set_xlim(left = 0)
ax.set_ylim(bottom = 0)
plt.title("Polycarb Force vs Displacement")
plt.ylabel('Load (N)')
plt.xlabel('Displacement (mm)')
plt.legend()


# Finding Sample Elastic slopes
def slopeFit(displacement, force, a, b):
    slope, intercept, R, P, Err = linregress(displacement[a:b], force[a:b])
    yP = 5 + next(i for i, x in enumerate(displacement) if force[i+5] < 0.95*slope*x + intercept)
    critLoad = force[yP]
    return slope, intercept, critLoad, R

for file in Files:
    if 'Ti 0deg' in file:
        Data[file]['slope'],Data[file]['intercept'],Data[file]['Critical Load (N)'],Data[file]['rval'] = slopeFit(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],0,100) #revise bounds later
    elif 'Ti 90deg' in file:
        Data[file]['slope'],Data[file]['intercept'],Data[file]['Critical Load (N)'],Data[file]['rval'] = slopeFit(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],0,100) #revise bounds later
    elif 'PC' in file:
        Data[file]['slope'],Data[file]['intercept'],test,Data[file]['rval'] = slopeFit(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],0,100) #revise bounds later
    else:
        Data[file]['slope'],Data[file]['intercept'],test,Data[file]['rval'] = slopeFit(Data[file]['Displacement (mm) '],Data[file]['Load (N)'],0,100) #revise bounds later


# Calculating max loads/critical load for plastic
for file in Files:
    Data[file]['Max Load (N)'] = max(Data[file]['Load (N)'])
    if file in np.concatenate([PC,PMMA]):
        Data[file]['Critical Load (N)'] = Data[file]['Max Load (N)']

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
    print(file[0:len(file)-4],'Critical Load (N):',Data[file]['Critical Load (N)'][0],'Conditional Crack Toughness:',Data[file]['Conditional Crack Toughness'][0],'Fracture Toughness:',Data[file]['Fracture Toughness'][0])

plt.show()















