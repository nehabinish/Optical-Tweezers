% cptrap.m - Counter-propagating trap in geometrical optics
%
% Optical forces arising from two counter-propagating rays impinging on a
% spherical particle within the geometrical optics approach. Only
% determinsitic optical forces are taken into account.
%
% See also RAY, PARTICLESPHERICAL, POINT, VECTOR.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters

% Particle and medium
R = 3e-6; % Particle radius [m]
np = 1.5; % Particle refractive index
nm = 1.33; % Medium refractive index
gamma = 6*pi*0.001*R; % Friction coefficient in water [N/s]
c = Point(0,-2e-6,-1e-6); % Initial position [m]

% Ray A
vA = Vector(-10e-6,0,0,4e-6,0,0); % Direction [m]
PA = 1e-3; % power [W]
polA = Vector(0,0,0,0,1,1); polA = vA*polA; polA = 1e-6*polA.versor(); % polarization
rA = Ray(vA,PA,polA);
        
% Ray B
vB = Vector(10e-6,0,0,-4e-6,0,0); % Direction [m]
PB = 1e-3; % power [W]
polB = Vector(0,0,0,0,1,1); polB = vB*polB; polB = 1e-6*polB.versor(); % polarization
rB = Ray(vB,PB,polB);

%% Simulation

for n = 1:1:100
    
    % Particle
    bead = ParticleSpherical(c,R,nm,np);

    % Scattered rays
    r_vecA = bead.scattering(rA,1e-18,0);
    r_vecB = bead.scattering(rB,1e-18,0);

    % Forces
    fA = bead.force(rA,1e-18,0);
    fB = bead.force(rB,1e-18,0);
    f = fA+fB;

    % Displacements
    Dt = 1e-3; % Timestep [S]
    DrA = fA*(Dt/gamma);
    DrB = fB*(Dt/gamma);
    Dr = DrA+DrB;
    
    % Center displacement
    c = Point(c.X+Dr.Vx,c.Y+Dr.Vy,c.Z+Dr.Vz);
    
    % Figure
    cla
    title(['T=' num2str(n*Dt) 's'])
    hold on
    bead.sp.plot('edgecolor','k');
    f = 1e+7*f; f.plot('color','r');
    rA.plot('color','g','linewidth',2);
    rB.plot('color','g','linewidth',2);
    r_vecA(1).t.plot('color','g','linewidth',2);
    r_vecB(1).t.plot('color','g','linewidth',2);
    r_vecA(2).t.plot('color','g','linewidth',2);
    r_vecB(2).t.plot('color','g','linewidth',2);
    hold off
    axis equal
    axis([-10e-6 10e-6 -5e-6 5e-6 -5e-6 5e-6])
    grid on
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    drawnow()

end