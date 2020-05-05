numOfSims = 10;
N = 1e7;
Dt = 1e-4;
sampleRate = 20;
R = 1e-6;
T  = 300;
eta = 0.001;
lx = 11e-6;
ly = 11e-6;
numOfParticles = 7;
wallShrink = 0;
baseFoldername = 'Parallel Harmonic ';
displayLive = false;

parfor i=1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    infoChamber(N, Dt, sampleRate, R, T, eta, lx, ly, numOfParticles, wallShrink, currFoldername, displayLive);
end

%% Combining the results from the different simulations into one heat map
load(strcat(baseFoldername,'1','/densityMatrix.mat'));
fullDensity = meanDensity;
for i = 2:1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    load(strcat(baseFoldername,num2str(i),'/densityMatrix.mat'));
    fullDensity = fullDensity + meanDensity;
end
fullDensity = fullDensity ./ numOfSims;
imagesc(fullDensity);
%% combining the results of the standard deviations to get a the characteristic correlation decay time
load(strcat(baseFoldername,'1','/densitySTDs.mat'));
fullSTDs = densitySTDs;
for i = 2:1:numOfSims
    currFoldername = strcat(baseFoldername,num2str(i));
    load(strcat(baseFoldername,num2str(i),'/densitySTDs.mat'));
    fullSTDs = fullSTDs + densitySTDs;
end
fullSTDs = fullSTDs ./ numOfSims;
t = 0:1/sampleRate:(length(fullSTDs)-1)*(1/sampleRate);
figure
plot(t,fullSTDs);
xlabel('time (sec)');
ylabel('Standard Deviation');
title('Averaged Standard Deviation of the particle densities');