classdef simConfig < handle
    properties
        numOfParticles % The amount of particles to simulate
        R % A (row) vector of the radii of the particles
        initPositions % The initial positions of the particles (2xnumOfParticles matrix)
        eta % The viscosity coefficient eta of the liquid
        T % The temperature at which the simulation is run
        N % The number of steps to run the simulation for
        Dt % The length (sec) of a time step to use in the simulation
        sampleRate % The rate (1/sec) at which the particle positions are recorded
        useWalls % Boolean determining whether to apply walls to the simulation
        wallPositionsX % Spatial position of the wall on the x axis
        wallPositionsY % Spacial position of the wall on the y axis
        useParticleRepulsion % Boolean determining whether particles apply hard core repulsion
        WCAEpsilon % The epsilon parameter for WCA repulsion
        useHydro % Determines whether to apply hydrodynamic interactions between the particles
        displayLive % A boolean to determine whether to show the particles live on screen during the simulation
    end
end