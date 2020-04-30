function [meanDensity, densitySTDs] = densityHeatMap(particlePositions, radii, xLimits, yLimits, pixelsPerLength, verbose, movie, saveFoldername)
% DENSITYHEATMAP  Generates a heat map of the densities of particles given
% their tracks
%   particlePositions:  The tracks of the particles used (numOfFrames,
%   numOfParticles, numOfDimensions)
%   radii: A vector containing the radii of the particles
%   xLimits, yLimits: The range in which to check the density (each is a
%   two-value vector of the lower and upper bounds)
%   pixelsPerLength: Determines how many pixels (matrix cells) would be
%   generated per length unit in the tracks (coarse-graining)
%   verbose: A boolean detemining whether to show the heat map being
%   generated live
%   movie: A boolean determining whether to save a movie of the evolution
%   of the heat map
%   saveFoldername: The folder in which to save the STD data, and the movie
%   if needed

    lx = abs(xLimits(2) - xLimits(1));
    ly = abs(yLimits(2) - yLimits(1));
    numOfFrames = size(particlePositions,1);
    numOfParticles = size(particlePositions, 2);
    
    % Creating a coarse-grained matrix of the tested area
    nx = round(pixelsPerLength.*lx);
    ny = round(pixelsPerLength.*ly);
    meanDensity = zeros(ny, nx);
    % Note that since these are images, we use the rows for the y axis and
    % the columns for the x axis, opposite from what the simulation uses.
    
    densitySTDs = zeros(numOfFrames, 1);
    radiiInMat = radii.*pixelsPerLength;
    figure(3);
    for frameInd = 1:numOfFrames
        for particleInd = 1:numOfParticles
            currPosition = squeeze(particlePositions(frameInd,particleInd,:));
            % Getting the center in of the particle in the coarse-grained
            % matrix
            positionInMat = [(currPosition(1) - xLimits(1)) .* pixelsPerLength,
                             (currPosition(2) - yLimits(1)) .* pixelsPerLength];
            % Finding the points in which the particle is present (a circle
            % of radius R around the center of the particle)
            cellsToRaise = getCirclePart(nx,ny,...
                                         positionInMat(1),positionInMat(2),...
                                         0, radiiInMat(particleInd),...
                                         0,2*pi);
            % Raising the intensity in the coarse grained matrix at the
            % points in which the particle is present.
            meanDensity(cellsToRaise) = meanDensity(cellsToRaise) + 1;
        end
        
        % Measuring the standard deviation of the mean intensity
        % distribition of all steps up to the current one. (The division by
        % the frame index is essentially normalization)
        densitySTDs(frameInd) = std(meanDensity(:) ./ frameInd);
        
        % Checking whether to show the data on screen
        if verbose || movie
            imagesc(meanDensity ./ frameInd);
            set(gca,'YDir','normal');
            title(frameInd);
            if movie
                % Note that F is not pre-allocated, since matlab forbids
                % pre-allocation of very large matrices
                F(frameInd) = getframe(gcf);
                drawnow
            end
        end
        pause(0.001);
    end
    % Saving the final matrix
    save(strcat(saveFoldername, '/densityMatrix.mat'),'meanDensity');
    % Saving the final figure
    saveas(gcf, strcat(saveFoldername,'/finalHeatMap.png'));
    % Saving the movie to the disk, if needed
    if movie
        saveMovie(F, 60, strcat(saveFoldername, '/heatMap.avi'));
    end
    % Saving the vector of density standard deviations
    save(strcat(saveFoldername, '/densitySTDs.mat'),'densitySTDs');
    % Presenting the plot of the standard deviation of the mean particle
    % density
    figure(4)
    plot(densitySTDs);
    xlabel('Number of steps');
    ylabel('Standard Deviation');
    title(strcat('StDev of the mean particle density ', saveFoldername));
    saveas(gcf, strcat(saveFoldername,'/STDs.png'));