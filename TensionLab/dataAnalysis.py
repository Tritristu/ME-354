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
    Area = Data[File]['Thickness(m)'][0]*Data[File]['Width(m)'][0]*(10**6)
    #Here's some example new dictionary calculations
    Data[File]['Strain (mm/mm)'] = Data[File]['Axial Strain (mm/mm)']  # this adds another strain column to the data
    Data[File]['Stress (MPa)'] = Data[File]['Load (N)']/Area #this adds stress to the data
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


#Youngs Modulus Calculation
#Here we will make a subplot to show a zoomed section
eZoom = 0.01; sZoom = 350
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
for File in Files:
    ax1.plot(Data[File]['Strain (mm/mm)'],Data[File]['Stress (MPa)'],label=File)
    ax2.plot(Data[File]['Strain (mm/mm)'],Data[File]['Stress (MPa)'],label=File)

#Plot the zoomed box
ax1.plot([0,0,eZoom,eZoom,0],[0,sZoom,sZoom,0,0],'r--')

#Add Labels
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax2.set_xlim(left = 0, right=eZoom)
ax2.set_ylim(bottom = 0, top=sZoom)
ax1.set_title("Engineering Stress vs Strain")
ax2.set_title("Engineering Stress vs Strain Magnified")
ax1.set_ylabel('Stress (MPa)')
ax2.set_ylabel('Stress (MPa)')
ax1.set_xlabel('Strain (mm/mm)')
ax2.set_xlabel('Strain (mm/mm)')
ax1.legend()
ax2.legend()


def modulusFit(Strain, Stress, a, b):
    '''This is a linear fit to data between the data indices for a and b. Note, this will
    return an error if a or b are outside the length of Strain and Stress.'''

    # Fit the modulus
    E, C, R, P, Err = linregress(Strain[a:b], Stress[
                                              a:b])  # The data outputs the slope (E), intercept (C), regression (R) value, P-value and standard error
    # Note: Python lets you save multivariable outputs with a comma, i.e. a,b=[1,2] will give a=1 and b=2

    # Make a line for the fit data
    Y = [0.0, max(Stress)]  # this is a list of length 2 for plotting the fit data later
    X = [(y - C) / E for y in
         Y]  # these are points that you can plot to visualize the data being fit, inverted from y=E*x+C, x=(y-C)/E
    return E, C, R, X, Y

#Checking Fit
fig = plt.figure(3)
ax = fig.gca()

# We'll test out two fit points and see the difference
a1 = 125
b1 = 200
a2 = 20
b2 = 150
for File in Files:
    # Save dummy variables to make the code cleaner below
    strain = Data[File]['Strain (mm/mm)'].values
    stress = Data[File]['Stress (MPa)'].values

    # Use the two fits
    E1, C1, R1, X1, Y1 = modulusFit(strain, stress, a1, b1)
    E2, C2, R2, X2, Y2 = modulusFit(strain, stress, a2, b2)

    # Plot the data
    ax.plot(strain, stress)
    ax.plot(strain[a1], stress[a1], 'rd')  # This is the first point we're fitting from
    ax.plot(strain[b1], stress[b1], 'rs')  # this is the last point we're fitting to
    #     ax.plot(strain[a1:b1],stress[a1:b1],'r.') #this will show all the data we're fitting
    #ax.plot(strain[a2], stress[a2], 'bd')
    #ax.plot(strain[b2], stress[b2], 'bs')
    #     ax.plot(strain[a2:b2],stress[a2:b2],'b.') #this will show all the data we're fitting

    # Plot the fits
    ax.plot(X1, Y1, label='Fit1, E=' + str(round(E1 * 1e-3, 1)) + ' GPa')
    #ax.plot(X2, Y2, label='Fit2, E=' + str(round(E2 * 1e-3, 1)) + ' GPa')

ax.set_xlim(left = 0, right=0.004)
ax.set_ylim(bottom = 0, top=400)
plt.title("Modulus Fit")
plt.ylabel('Stress (MPa)')
plt.xlabel('Strain (mm/mm)')
plt.legend()


#FIND WHAT VALUES OF A AND B ARE USED BY OTHER PEOPLE OR FIND THEM YOURSELF
a = [125,125,125,125]
b = [200,200,200,200]
youngsModuli = []
avgYoungsModuli = []
stdYoungsModuli = []
ER2Values = []
inc = 0
inc2 = 0
for File in Files:

    # Save dummy variables to make the code cleaner below

    strain = Data[File]['Strain (mm/mm)'].values
    stress = Data[File]['Stress (MPa)'].values

    # Use the two fits
    E,C,R,X,Y = modulusFit(strain,stress,a[int(inc/3)],b[int(inc/3)])
    youngsModuli += [E]
    #print(File)
    ER2Values += [R**2]
    inc = inc+1
for i in range(4):
    avgYoungsModuli+= [mean(youngsModuli[i*3:(i+1)*3])/1000]
    stdYoungsModuli+= [np.std(youngsModuli[i*3:(i+1)*3])/1000]

print(avgYoungsModuli)
print(stdYoungsModuli)

#finding poisons ratio

def PoissonFit(axialStrain, transverseStrain, a, b):
    '''This is a linear fit to data between the data indices for a and b. Note, this will
    return an error if a or b are outside the length of Strain.'''

    # Fit the modulus
    nu, C, R, P, Err = linregress(axialStrain[a:b], transverseStrain[
                                                    a:b])  # The data outputs the slope (nu), intercept (C), regression (R) value, P-value and standard error
    # Note: Python lets you save multivariable outputs with a comma, i.e. a,b=[1,2] will give a=1 and b=2

    # Make a line for the fit data
    Y = [0.0, transverseStrain[round(1.5 * b)]]
    X = [(y - C) / nu for y in
         Y]  # these are points that you can plot to visualize the data being fit, inverted from y=nu*x+C, x=(y-C)/nu
    return nu, R, X, Y


fig = plt.figure(5)
ax = fig.gca()
a = 0;
b = 10  # Note: you will have to play with these values for a given test type
for File in Files:
    # Create dummy variables to make plotting easier
    aStrain = Data[File]['Axial Strain (mm/mm)'].values
    tStrain = Data[File]['Transverse Strain (mm/mm)'].values

    # Do the fit
    nu, R, X, Y = PoissonFit(aStrain, tStrain, a, b)

    # Plot things
    ax.plot(aStrain, tStrain)
    ax.plot(aStrain[a], tStrain[a], 'rd')  # this is the first point we're fitting from
    ax.plot(aStrain[b], tStrain[b], 'bs')  # this is the last point we're fitting from
    ax.plot(aStrain[a:b], tStrain[a:b], 'b.')  # These are all the points we're fitting

    # Plot the fits
    ax.plot(X, Y, label=r"Poisson's Ratio Fit, $\nu$=" + str(round(nu, 3)))

# Format the plot
ax.set_xlim(left=0, right=0.01)
ax.set_ylim(bottom=0, top=0.002)
plt.title("Axial vs Transverse Strain Magnified")
plt.xlabel('Axial Strain (mm/mm)')
plt.ylabel('Transverse Strain (mm/mm)')
plt.legend()

#FIND WHAT VALUES OF A AND B ARE USED BY OTHER PEOPLE OR FIND THEM YOURSELF


a = [10,20,20,10]
b = [100,200,250,350]
PoissonsRatios = []
nuR2Values = []
avgPoissonsRatio = []
avgnuR2Values = []
inc = 0
for File in Files:
    # Save dummy variables to make the code cleaner below
    aStrain = Data[File]['Axial Strain (mm/mm)'].values
    tStrain = Data[File]['Transverse Strain (mm/mm)'].values

    # Use the two fits
    nu,R,X,Y = PoissonFit(aStrain,tStrain,a[int(inc/3)],b[int(inc/3)])

    PoissonsRatios += [nu]
    nuR2Values += [R**2]
    inc = inc+1
for i in range(4):
    avgPoissonsRatio += [mean(PoissonsRatios[i * 3:(i + 1) * 3])]
    avgnuR2Values += [std(nuR2Values[i * 3:(i + 1) * 3])]

plt.show()

