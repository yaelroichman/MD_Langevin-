function f=getLJForce(r, epsilon, sigma)
%Calculate force from Lennard-Jones potential
f = 48*epsilon*((1/r)*((sigma/r)^12-0.5*(sigma/r)^6));