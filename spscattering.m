% spscattering.m - Scattering of a ray on a spherical particle
%
% Calculation of the scattering efficiencies corresponding to the
% scattering of a ray on a spherical particle as a function of the
% incidence angle.
%
% See also RAY, SPHERICALPARTICLE.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters

% Particle and medium
R = 1e-6;  % Particle radius [m]
np = 1.50;  % Particle refractive index
nm = 1.33;  % Medium refractive index

%% Simulation

% Particle
c = Point(0,0,0);  % Particle center [m]
bead = ParticleSpherical(c,R,nm,np);

% Rays
theta = [0:1:89.9]/180*pi;  % Incidence angles [rad]
v = Vector(-2*R*ones(size(theta)),zeros(size(theta)),R*sin(theta),ones(size(theta)),zeros(size(theta)),zeros(size(theta))); % Direction
P = ones(size(theta));  % Power [W]
pol = Vector(zeros(size(theta)),zeros(size(theta)),zeros(size(theta)),zeros(size(theta)),ones(size(theta)),ones(size(theta))); pol = v*pol; pol = pol.versor(); % Polarization
r = Ray(v,P,pol);

% Scattering coefficients
f = bead.force(r,1e-18,100);
cm = PhysConst.c0/bead.nm;  % Speed of light in medium [m/s]
Q = f./(r.P/cm);  % Scattering efficiency

% Figure
figure

title('Exact Trapping Efficiencies')
hold on
plot(theta/pi*180,Q.Vx,'r-.','linewidth',2.5)
plot(theta/pi*180,Q.Vz,'b--','linewidth',2.5)
plot(theta/pi*180,sqrt(Q.Vx.^2+Q.Vz.^2),'k-','linewidth',2.5)
hold off
legend('Q_s','Q_g','Q','Location','NorthWest')
box on
grid on
xlabel('Incidence angle')
ylabel('Trapping efficiencies')