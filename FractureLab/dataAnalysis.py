import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # use seaborne instead?
from os import listdir
from scipy.stats import linregress

# Constants


# reading data files
Files = [x for x in listdir('FractureLab') if '.csv' in x]
Ti0 = [x for x in listdir('FractureLab') if 'Ti 0deg' in x]
Ti90 = [x for x in listdir('FractureLab') if 'Ti 90deg' in x]
PMMA = [x for x in listdir('FractureLab') if 'PMMA' in x]
PC = [x for x in listdir('FractureLab') if 'PC' in x]

Data = {x:{} for x in Files}
for File in Files:
    Data[File] = pd.read_csv('FractureLab/' + File)

# for file in [PC,PMMA]:
    # Data[file]['Max Load (N)'] = max(Data[file]['Load (N)'])


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

# plt.show()















