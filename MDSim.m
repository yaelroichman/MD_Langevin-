function particlePositions = MDSim(cfg, forcesFunc, feedbackCheckFunc, feedbackFunc, printFunc, addedData)
%% Basic definitions
close all
kB = physconst('Boltzmann'); % Boltzmann constant [J/K]
R = cfg.R(1); % Currently works only with one constant R
gamma = 6*pi*R*cfg.eta; % friction coefficient
D = kB*cfg.T/gamma; % diffusion coefficient
d = 2; % The dimension of the problem. Currently ONLY WORKS FOR 2d!!
samplePeriod = round(1 / (cfg.Dt*cfg.sampleRate)); % Defines the samplePeriod at which we will sample the simulation. 
particlePositions = zeros(1+ceil(cfg.savePeriod/samplePeriod),cfg.numOfParticles, d); % Preallocate space
particlePositions(1,:,:) = cfg.initPositions'; % sets the first step as defined in configuration
%% Checking if the save directory already exists
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
%% Saving the configuration
save(strcat(cfg.saveFoldername, '/cfg.mat'), 'cfg');
%% Setting up the run variables
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

%% Plotting the initial placements
printFunc(cfg, currStepData, addedData);
xlabel('x [m]');
ylabel('y [m]');
title('Initial placement');
pause(0.01);
%% Checking whether to run hydrodynamic interactions
if ~cfg.useHydro
    Dx = ones(cfg.numOfParticles,1).*D;
    Dy = Dx;
    Ax = ones(cfg.numOfParticles,1).*sqrt(2*D);
    Ay = Ax;
end
%% Running the simulation
sampleInd = 2;
for i = 2:1:cfg.N
    currStepData.stepNum = i;
    %% Printing data
    if mod(i, samplePeriod) == 0
        if cfg.displayLive
            printFunc(cfg, currStepData, addedData);
        end
        
        particlePositions(sampleInd,:,:) = currStepData.particlePositions;
        sampleInd = sampleInd + 1;
    end
    %% Saving the steps according to the save period parameter
    if mod(i, cfg.savePeriod) == 0 || i == cfg.N
        %% Saving the latest particle positions in .csv files
        dlmwrite(strcat(cfg.saveFoldername, '/pos_x.csv'),...
            particlePositions(1:sampleInd-1,:,1),...
            '-append');
        dlmwrite(strcat(cfg.saveFoldername, '/pos_y.csv'),...
            particlePositions(1:sampleInd-1,:,2),...
            '-append');
        
        particlePositions(1:sampleInd-1,:,:) = 0;

        %% Saving the additional tracked data in a .mat file.
        if exist(strcat(cfg.saveFoldername, '/data.mat'), 'file')
            delete(strcat(cfg.saveFoldername, '/data.mat'));
        end
        save(strcat(cfg.saveFoldername, '/data.mat'), 'addedData');
        
        sampleInd = 1;
        100*(i/cfg.N)
    end
    %% Forces computation
    [fx, fy] = ...
        forcesFunc(cfg, currStepData, addedData);
    %% Hydrodynamic interactions computation
    if cfg.useHydro
        DMat = rotnePrager(currStepData.particlePositions,R, D);
        rootMat = chol(DMat,'lower');
        DVec = DMat*ones(cfg.numOfParticles*d,1);
        AVec = rootMat*ones(cfg.numOfParticles*d,1);
        Dx = DVec(1:2:end);
        Dy = DVec(2:2:end);
        Ax = AVec(1:2:end);
        Ay = AVec(2:2:end);
    end
    %% Running the step
    currStepData.particlePositions(:,1) = currStepData.particlePositions(:,1) +...
        Ax.*sqrt(cfg.Dt).*randn(cfg.numOfParticles,1) +...
        (Dx./(kB.*cfg.T)).*fx.*cfg.Dt;
    currStepData.particlePositions(:,2) = currStepData.particlePositions(:,2) +...
        Ay.*sqrt(cfg.Dt).*randn(cfg.numOfParticles,1) +...
        (Dy./(kB.*cfg.T)).*fy.*cfg.Dt;
    %% Checking whether to apply a feedback to the system, and applying it if necessary
    if feedbackCheckFunc(cfg, currStepData, addedData)
        feedbackFunc(cfg, currStepData, addedData);
    end
end

%% Reading the output file
particlePositionsX = dlmread(strcat(cfg.saveFoldername, '/pos_x.csv'));
particlePositionsY = dlmread(strcat(cfg.saveFoldername, '/pos_y.csv'));
particlePositions = zeros(size(particlePositionsX, 1), cfg.numOfParticles, 2);
particlePositions(:,:,2) = particlePositionsY;
particlePositions(:,:,1) = particlePositionsX;