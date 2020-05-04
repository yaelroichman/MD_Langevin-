classdef simConfig < handle
    properties
        xlimits %set borders in field of view- x axis
        ylimits %set borders in field of view- y axis
        numOfParticles % The amount of particles to simulate
        R % A (row) vector of the radii of the particles
        initPositions % The initial positions of the particles (2xnumOfParticles matrix)
        eta % The viscosity coefficient eta of the liquid
        T % The temperature at which the simulation is run
        N % The number of steps to run the simulation for
        Dt % The length (sec) of a time step to use in the simulation
        sampleRate % The rate (1/sec) at which the particle positions are recorded
        
        useWalls % Boolean determining whether to apply walls to the simulation
        wallRepulsionType % A string determining what kind of repulsion the walls apply (WCA/Harmonic/Gaussian)
        wallHarmonicK % The repulsion constant K for harmonic wall repulsion (Joule/meter)
        wallGaussianA % The amplitude for Gaussian wall repulsion
        wallGaussianS % The STD for gaussian wall repulsion
        wallPositionsX % Spatial position of the wall on the x axis
        wallPositionsY % Spacial position of the wall on the y axis
        
        useParticleRepulsion % Boolean determining whether particles apply hard core repulsion
        WCAEpsilon % The epsilon parameter for WCA repulsion
        
        useTraps % Boolean determining whether to use optical traps
        initTrapPositions % The initial positions of the particles (2xnumOfTrap matrix)
        A % Trap depth can be a scalar (if all traps are the same) or a (2xnumOfTrap matrix) 
        s % Trap size can be a scalar (if all traps are the same) or a (2xnumOfTrap matrix) 
        
        useHydro % Determines whether to apply hydrodynamic interactions between the particles
        displayLive % A boolean to determine whether to show the particles live on screen during the simulation
        saveFoldername % The name of the file to save the data into
        savePeriod % The amount of steps between writes to the disk
    end
end