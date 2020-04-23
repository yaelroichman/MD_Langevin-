
# MD_Langevin
### Matlab based Stokesian and Brownian dynamics code for simulating colloidal motion in fluid

## Purpose

To calculate the trajectories of colloidal particles diffusing in a fluid. 

The code can:
1. Take into account hydrodynamic interactions between the particles under the stokes approximation with the Rotne–Prager correction. 
2. Take into account body forces and interaction forces (e.g. hard core repulsion) acting on the diffusing particles.


## MDSim
### The main body of the simulation. Will be explained in detail followed by an explanation about the function that go into the sim. 
```matlab
particlePositions = MDSim(cfg, forcesFunc, feedbackCheckFunc, feedbackFunc, printFunc, addedData)
```
#### Dissecting MDSim
All, non standard functions, can be found in the ["Function explanation"](#function-explanation) section.

##### %% Basic definitions
We start from basic defintions that pull from the configuration file (cfg.something).
When preallocating space we add 1 for the initial positions, and we use ceil just in case the `samplePeriod` isn't a factor of the `cfg.savePeriod`.
```matlab
close all
kB = physconst('Boltzmann'); % Boltzmann constant [J/K]
R = cfg.R(1); % Currently works only with one constant R
gamma = 6*pi*R*cfg.eta; % friction coefficient
D = kB*cfg.T/gamma; % diffusion coefficient
d = 2; % The dimension of the problem. Currently ONLY WORKS FOR 2d!!
samplePeriod = round(1 / (cfg.Dt*cfg.sampleRate)); % Defines the samplePeriod at which we will sample the simulation. 
particlePositions = zeros(1+ceil(cfg.savePeriod/samplePeriod),cfg.numOfParticles, d); % Preallocate space
particlePositions(1,:,:) = cfg.initPositions'; % sets the first step as defined in configuration
```
***Notice:*** see that it only works in 2d currently.
 
The `particlePositions(1,:,:)` changes our `cfg.initPositions'` from being a `[2XNumOfParticles]` into lines of ***X*** positions and ***Y*** positions on different 3rd dimension. 
##### Example 1
from `initPos = [0 0 ; 1 1 ; 2 2];`
*initPos =*

     0     0
     1     1
     2     2`
to `partPos(1,:,:) = initPos';`
*ans(:,:,1) =*

     0     1     2


*ans(:,:,2) =*

     0     1     2

---

##### %% Checking if the save directory already exists
Creates / Cleans the folder where we will save.
```matlab
if ~exist(cfg.saveFoldername, 'dir')
    mkdir(cfg.saveFoldername);
else
    answer = questdlg('Save folder already exists. Overwrite?');
    if strcmp(answer, 'Yes')
        delete(strcat(cfg.saveFoldername, '/pos_x.csv'));
        delete(strcat(cfg.saveFoldername, '/pos_y.csv'));
        delete(strcat(cfg.saveFoldername, '/cfg.m'));
    else
        error('Aborting sim to avoid overwriting the existing saved file');
    end
end
```
##### %% Saving the configuration
Saves the configuration file into the save folder.
```matlab
save(strcat(cfg.saveFoldername, '/cfg.m'), 'cfg');
```
---

##### %% Setting up the run variables
`simStepData` is like out configuration file, but only for things that may change every step (e.g. particlePositions, step number, wall positions (*if we chose to use walls*)).
You will always have the `particlePositions` and the `stepNum` (defined later), available at each simulation step.
```matlab
currStepData = simStepData;
currStepData.particlePositions = (squeeze(particlePositions(1,:,:)));
if cfg.useWalls
    currStepData.wallPositionsX = cfg.wallPositionsX;
    currStepData.wallPositionsY = cfg.wallPositionsY;
end

if cfg.useTraps
    currStepData.A = cfg.A;
    currStepData.s = cfg.s;
    trapPositions(1,:,:) = cfg.initTrapPositions';
    currStepData.trapPositions = (squeeze(trapPositions(1,:,:)));
end
```
***Notice:*** Currently the `simStepData` contains the option to add wall positions, and trap position information. It is added here because our simulations will probably need to move walls or change trap positions / intensity during the run of the simulation. 
When adding new forces or feedback functions you most likely will not need to add things to the `simStepData`.
*In any case you can refer to the* **["How to personalize"](#how-to-personalize)** *section for any such concerns.*

---

##### %% Plotting the initial placements
prints the initial positions as defined in configuration file.
```matlab
printFunc(cfg, currStepData, addedData);
xlabel('x [m]');
ylabel('y [m]');
title('Initial placement');
pause(0.01);
```

---

##### %% Checking whether to run hydrodynamic interactions
If we don't want to use hydrodynamic interactions there is no need to calculate them during the simulation. We define our Diffusion vector and its square root appropriately to fit within our simulation.
D is taken from  *[%% Basic definitions](#basic-definitions)*.
```matlab
if ~cfg.useHydro
    Dx = ones(1, numOfParticles).*D;
    Dy = Dx;
    Ax = ones(1, numOfParticles).*sqrt(2*D);
    Ay = Ax;
end
```

---

#### %% Running the simulation
Runs from the second place till the defined `cfg.N` steps. Also defines the step number.
```matlab
sampleInd = 2;
for i = 2:1:cfg.N
    currStepData.stepNum = i
```

##### %% Printing data
Saves the current `particlePositions`. If `cfg.displayLive` is True will go to the defined `printFunc`.
```matlab
    if mod(i, samplePeriod) == 0
        if cfg.displayLive
            printFunc(cfg, currStepData, addedData);
        end
        
        particlePositions(sampleInd,:,:) = currStepData.particlePositions;
        sampleInd = sampleInd + 1;
    end
```

##### %% Saving the steps according to the save period parameter
When the simulation reaches the end of a `savePeriod` or the end of the defined number of steps (`cfg.N`) it appends the positions to the relevant file.
Resets the `ParticlePositions` matrix, and the `SampleInd`.
## Check with gilad
outputs completion precentage.
```matlab
    if mod(i, cfg.savePeriod) == 0 || i == cfg.N
        dlmwrite(strcat(cfg.saveFoldername, '/pos_x.csv'),...
            particlePositions(1:sampleInd-1,:,1),...
            '-append');
        dlmwrite(strcat(cfg.saveFoldername, '/pos_y.csv'),...
            particlePositions(1:sampleInd-1,:,2),...
            '-append');
            
        particlePositions(1:sampleInd-1,:,:) = 0;
        
        if exist(strcat(cfg.saveFoldername, '/data.m'), 'file')
            delete(strcat(cfg.saveFoldername, '/data.m'));
        end
        save(strcat(cfg.saveFoldername, '/data.m'), 'addedData');
        
        sampleInd = 1;
        100*(i/cfg.N)
    end
```

##### %% Forces computation
Uses the gives `forcesFunc` to compute the forces in both axes.
If we chose to use hydrodynamic interactions uses `rotnePrager` to compute the Diffusion matrix [`DMat`], then uses Cholesky decomposition to get something like the "square root" of the matrix [`rootMat`] for easier calculations after [[1-3](#references)]. Turns the matrices into vectors of different axes which will be used later.
```matlab    
    [fx, fy] = ...
        forcesFunc(cfg, currStepData, addedData);
    %% Hydrodynamic interactions computation
    if cfg.useHydro
        DMat = rotnePrager(currStepData.particlePositions,R, D);
        rootMat = chol(DMat,'lower');
        DVec = DMat*ones(cfg.numOfParticles.*d,1);
        AVec = rootMat*ones(cfg.numOfParticles.*d,1);
        Dx = DVec(1:2:end);
        Dy = DVec(2:2:end);
        Ax = AVec(1:2:end);
        Ay = AVec(2:2:end);
    end
```

##### %% Running the step
Calculates the current position of each particle (on each axis) using a random white noise, and the forces from the `forcesFunc` [[4](#references)].
```matlab 
    currStepData.particlePositions(:,1) = currStepData.particlePositions(:,1) +...
        Ax.*sqrt(cfg.Dt).*randn(cfg.numOfParticles,1) +...
        (Dx./(kB.*cfg.T)).*fx.*cfg.Dt;
    currStepData.particlePositions(:,2) = currStepData.particlePositions(:,2) +...
        Ay.*sqrt(cfg.Dt).*randn(cfg.numOfParticles,1) +...
        (Dy./(kB.*cfg.T)).*fy.*cfg.Dt;
```

##### %% Checking whether to apply a feedback to the system, and applying it if necessary
Checks against some `feedbackCheckFunc` given by the user if to use some `feedbackFunc` also given by the user.
```matlab
    if feedbackCheckFunc(cfg, currStepData, addedData)
        feedbackFunc(cfg, currStepData, addedData);
    end
end
```

##### %% Reading the output file
After the simulation finishes it reads back into matlab the full `particlePositions` matrice.
```matlab
particlePositionsX = dlmread(strcat(cfg.saveFoldername, '/pos_x.csv'));
particlePositionsY = dlmread(strcat(cfg.saveFoldername, '/pos_y.csv'));
particlePositions = zeros(size(particlePositionsX, 1), cfg.numOfParticles, 2);
particlePositions(:,:,2) = particlePositionsY;
particlePositions(:,:,1) = particlePositionsX;
```

---
## Function explanation
cfg - input
simStepdata - calls different sections between MDSim
addedData - optional tank of data

forcesFunc

rotnePrager

 feedbackCheckFunc
 feedbackFunc
 printFunc
 

## Basic examples

Please read [Examples.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for some basic usage examples

## How to personalize

Please read [Personalization.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details about how to add new forces or other functions.

## References

1. J. F. Brady and G. Bossis, Ann. Rev. Fluid Mech. 20, 111 (1988).
2. Y. Sokolov, D. Frydel, D. G. Grier, H. Diamant, and Y. Roichman, Phys. Rev. Lett., 107, 158302, (2011).
3. H. Nagar and Y. Roichman, Phys. Rev. E. 90, 042302 (2014).
4. G. Volpe and G. Volpe, American Journal of Physics 81, 224 (2013);
