% otgo.m - Optical tweezers in geometrical optics
%
% Simulation of an optical tweezers using geometrical optics.
% The optical force on a spherical particle is calculated as a function of
% the particle postion with respect to the focal point.
%
% See also BEAMGAUSS, RAY, SPHERICALPARTICLE, POINT, VECTOR.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialisation of the workspace
clear all;
close all;
clc;

%% Parameters

% Particle and medium
R = 3e-6;       % Particle radius [m]
np = 1.5;       % Particle refractive index
nm = 1.33;      % Medium refractive index

% Focusing
f = 100e-6;     % Focal length [m]
NA = 1.3;       % Numerical aperture
L = f*NA/nm;    % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4;     % x electric field [V/m]
Ey0 = 1i*1e+4;  % y electric field [V/m]
w0 = 5e-3;      % Beam waist [m]
Nphi = 16;      % Azimuthal divisions
Nr = 16;        % Radial divisions
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);

%% Simulation

% Calculates set of rays corresponding to optical beam
r = Ray.beam2focused(bg,f);

% Executes calculation of force as a function of particle position
for x = [-3e-6:.5e-6:3e-6]  % x position [m]
    for y = [-3e-6:.5e-6:3e-6]  % y position [m]
        for z = [-6e-6:.5e-6:6e-6]  % z position [m]
                
            % Particle
            bead = ParticleSpherical(Point(x,y,z),R,nm,np);
            
            % Forces
            forces = bead.force(r);
            force = Vector(x,y,z, ...
                sum(forces.Vx(isfinite(forces.Vx))), ...
                sum(forces.Vy(isfinite(forces.Vy))), ...
                sum(forces.Vz(isfinite(forces.Vz))) ...
                );

            % Figure
            cla
            title(['x=' num2str(x*1e+6) 'um y=' num2str(y*1e+6) 'um z=' num2str(z*1e+6) 'um'])
            hold on
            plot3(0,0,0,'.r');
            r.versor.plot('color','g','scale',[1 1e-5]);
            bead.sp.plot('facecolor','k','edgecolor','k');
            force = 1e+7*force; 
            force.plot('color','r');
            hold off
            axis equal
            grid on
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            drawnow()

        end
    end
end