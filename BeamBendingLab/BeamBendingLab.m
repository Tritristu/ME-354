clear;close all;clc;

% 6061-T6 Aluminum Material Constants
ELASTIC_MODULUS = 69E9;
yieldStrength = 275E6;
poissonsRatio = 0.33;

% 3-Point Beam Dimensional Data
LENGTH = 700E-3;
WIDTH_3PT = 38.38E-3;
HEIGHT_3PT = 25.75E-3;
MOMENT_OF_INERTIA_3PT = (WIDTH_3PT * HEIGHT_3PT.^3)/12;

Z_POSITIONS_3PT = [-HEIGHT_3PT/2;
                    (12.25E-3)-(HEIGHT_3PT/2);
                    HEIGHT_3PT/2;
                   -HEIGHT_3PT/2]';

% Moment is associated with x position, so I included a scale for each
% gauge (bad practice, should have done them in different steps but am layzy)
WEIGHTED_Z_POSITIONS_3PT = [(-HEIGHT_3PT/2).*(LENGTH./2);
               ((12.25E-3)-(HEIGHT_3PT/2)).*(LENGTH-0.380);
               (HEIGHT_3PT/2).*(LENGTH./2);
               (-HEIGHT_3PT/2).*(LENGTH-0.450)]';

% 4-Point Beam Dimensional Data
WIDTH_4PT = 38.36E-3;
HEIGHT_4PT = 25.68E-3;
OFFSET = 200E-3;
MOMENT_OF_INERTIA_4PT = (WIDTH_4PT * HEIGHT_4PT.^3)/12;
Z_POSITIONS_4PT = [-HEIGHT_4PT/2;
                            (12.25E-3)-(HEIGHT_4PT/2);
                            HEIGHT_4PT/2;
                            -HEIGHT_4PT/2]';

% Data
realXStrains3Pt = [-224,4,238,-159;
                   -456,1,471,-326;
                   -688,3,704,-492;
                   -918,-5,936,-657].*1E-6;
realYStrains3Pt = [95,9,-79, NaN;
                   187,10,-164,NaN;
                   277,10,-248,NaN;
                   367,13,-332,NaN].*1E-6;

real45Rotations3Pt = [-67,7,90, -49;
                      -139,9,171,-106;
                      -208,12,250,-163;
                      -274,19,325,-212].*1E-6;

realXStrains4Pt = [-133,1,135,-141;
                   -264,1,268,-280;
                   -395,1,401,-419;
                   -528,1,535,-559].*1E-6;
realYStrains4Pt = [50,0,-48, NaN;
                   97,1,-95,NaN;
                   146,2,-142,NaN;
                   194,2,-190,NaN].*1E-6;

real45Rotations4Pt = [-30,11,31, -62;
                      -66,17,68,-118;
                      -105,19,110,-170;
                      -145,21,151,-221].*1E-6;


loads = [400;800;1200;1600];
actualDisplacements3Pt = [0.99;1.93;2.84;3.71];
actualDisplacements4Pt = [0.74;1.42;2.08;2.74];

shearForces3Pt = loads./2;
shearStresses3Pt = (shearForces3Pt./(2.*MOMENT_OF_INERTIA_3PT)).*((HEIGHT_3PT^2/4)-(Z_POSITIONS_3PT).^2);
shearStressesU3Pt = 0.79.*shearStresses3Pt;
shearStrains3Pt = 2.*shearStresses3Pt.*(1+poissonsRatio)./ELASTIC_MODULUS;
shearStrainsU3Pt = 0.79.*shearStrains3Pt;


Moments3Pt = loads./2;
xStresses3Pt = (Moments3Pt.*WEIGHTED_Z_POSITIONS_3PT)./MOMENT_OF_INERTIA_3PT;
xStressesU3Pt = 0.31.*xStresses3Pt;
xStrains3Pt = xStresses3Pt./ELASTIC_MODULUS;
yStrains3Pt = -(xStresses3Pt.*poissonsRatio)./ELASTIC_MODULUS;
xStrainsU3Pt = 0.31.*xStrains3Pt;
yStrainsU3Pt = 0.31.*yStrains3Pt;


Moments4Pt = [400/2;800/2;1200/2;1600/2].*OFFSET;
xStresses4Pt = (Moments4Pt.*Z_POSITIONS_4PT)./MOMENT_OF_INERTIA_4PT;
xStressesU4Pt = 0.25.*xStresses4Pt;
xStrains4Pt = xStresses4Pt./ELASTIC_MODULUS;
yStrains4Pt = -(xStresses4Pt.*poissonsRatio)./ELASTIC_MODULUS;
xStrainsU4Pt = 0.25.*xStrains4Pt;
yStrainsU4Pt = 0.25.*yStrains4Pt;


MaxDisplacement3Pt = @(load) (load.*LENGTH.^3)./(48.*ELASTIC_MODULUS.*MOMENT_OF_INERTIA_3PT);
MaxDisplacement4Pt = @(load) ((load.*OFFSET).*(3*LENGTH^2-4*OFFSET^2))./(48*ELASTIC_MODULUS*MOMENT_OF_INERTIA_4PT);

displacements3Pt = MaxDisplacement3Pt(loads);
displacements4Pt = MaxDisplacement4Pt(loads);

strainRotationEq = @(xStrain,Ystrain,shear,angle) ((xStrain + Ystrain)./2)+((xStrain - Ystrain)./2).*cosd(2.*angle) + (shear./2).*sind(2.*angle);
rotatedStrains3Pt = strainRotationEq(xStrains3Pt,yStrains3Pt,shearStrains3Pt,45);
rotatedStrains4Pt = strainRotationEq(xStrains4Pt,yStrains4Pt,0,45);


errorCalc = @(theoretical, real) abs((real-theoretical)./theoretical).*100;
XStrainError3Pt = errorCalc(xStrains3Pt,realXStrains3Pt);
YStrainError3Pt = errorCalc(yStrains3Pt,realYStrains3Pt);
rotatedError3Pt = errorCalc(rotatedStrains3Pt,real45Rotations3Pt);
reverseRotatedError3Pt = errorCalc(rotatedStrains3Pt(:,4),[-41;-92;-143;-203].*1E-6);

XStrainError4Pt = errorCalc(xStrains4Pt,realXStrains4Pt);
YStrainError4Pt = errorCalc(yStrains4Pt,realYStrains4Pt);
rotatedError4Pt = errorCalc(rotatedStrains4Pt,real45Rotations4Pt);
reverseRotatedError4Pt = errorCalc(rotatedStrains3Pt(:,4),[-30;-67;-107;-148].*1E-6);

%percentErrorU =@(theoretical,real,theoCert,realCert) sqrt((sqrt((realCert).^2 + (theoCert).^2)./(theoretical-real)).^2 +(theoCert./theoretical).^2);
%XStrainErrorU3Pt = percentErrorU(xStrains3Pt,realXStrains3Pt,xStrainsU3Pt,ones(4,4).*0.5E-6);
%YStrainErrorU3Pt = percentErrorU(yStrains3Pt,realYStrains3Pt,yStrainsU3Pt,ones(4,4).*0.5E-6);
%rotatedErrorU3Pt = ercentErrorU(rotatedStrains3Pt,real45Rotations3Pt,yStrainsU3Pt,ones(4,4).*0.5E-6);



%% Graphing

figure(1)
plot(loads,displacements3Pt.*1000, '-b',LineWidth=2);
hold on;
plot(loads,actualDisplacements3Pt, 'hexagramr',LineWidth=2);
hold on;
errorbar(loads, actualDisplacements3Pt, 0.05.*ones(1,4),'r.');
title("Three Point Bending Midpoint Displacements");
xlabel('Applied Load (N)');
ylabel('Displacement (mm)');
xlim([300,1700]);
legend('Theoretical Displacements', 'Real Displacements','Location','northwest');


figure(2)
plot(loads,displacements4Pt.*1000, '-b',LineWidth=2);
hold on;
plot(loads,actualDisplacements4Pt, 'hexagramr',LineWidth=2);
hold on;
errorbar(loads, actualDisplacements4Pt, 0.05.*ones(1,4),'r.');
title("Four Point Bending Midpoint Displacements");
xlabel('Applied Load (N)');
ylabel('Displacement (mm)');
xlim([300,1700]);
legend('Theoretical Displacements', 'Real Displacements','Location','northwest');

figure(3)

subplot(2,2,1)
plot(loads,XStrainError3Pt(:,1), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError3Pt(:,1), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError3Pt(:,1), 'gv',LineWidth=1.5);
subtitle("\bf Strain Rosette 1");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,30]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,2)
plot(loads,XStrainError3Pt(:,2), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError3Pt(:,2), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError3Pt(:,2), 'gv',LineWidth=1.5);
subtitle("\bf Strain Rosette 2");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,250]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,3)
plot(loads,XStrainError3Pt(:,3), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError3Pt(:,3), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError3Pt(:,3), 'gv',LineWidth=1.5);
subtitle("\bf Strain Rosette 3");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,15]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,4)
plot(loads,reverseRotatedError3Pt, 'bo',LineWidth=1.5);
hold on;
plot(loads,XStrainError3Pt(:,4), 'rx',LineWidth=1.5);
hold on;
plot(loads,rotatedError3Pt(:,4), 'gv',LineWidth=1.5);
subtitle("\bf Strain Rosette 4");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,40]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

sgtitle('\bf 3-Point Beam Bending Errors');


figure(4)

subplot(2,2,1)
plot(loads,XStrainError4Pt(:,1), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError4Pt(:,1), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError4Pt(:,1), 'gv',LineWidth=1.5);

subtitle("\bf Strain Rosette 1");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,40]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,2)
plot(loads,XStrainError4Pt(:,2), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError4Pt(:,2), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError4Pt(:,2), 'gv',LineWidth=1.5);

subtitle("\bf Strain Rosette 2");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,700]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,3)
plot(loads,XStrainError4Pt(:,3), 'bo',LineWidth=1.5);
hold on;
plot(loads,rotatedError4Pt(:,3), 'rx',LineWidth=1.5);
hold on;
plot(loads,YStrainError4Pt(:,3), 'gv',LineWidth=1.5);
subtitle("\bf Strain Rosette 3");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,35]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

subplot(2,2,4)
plot(loads,reverseRotatedError4Pt, 'bo',LineWidth=1.5);
hold on;
plot(loads,XStrainError4Pt(:,4), 'rx',LineWidth=1.5);
hold on;
plot(loads,rotatedError4Pt(:,4), 'gv',LineWidth=1.5);

subtitle("\bf Strain Rosette 4");
xlabel('\it Applied Load (N)');
ylabel('\it Error (%)');
xlim([300,1700]);
ylim([0,55]);
legend('\epsilon_{a} errors', '\epsilon_{b} errors','\epsilon_{c} errors','Location','northeast');

sgtitle('\bf 4-Point Beam Bending Errors');

