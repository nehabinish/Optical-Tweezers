% focus.m - Focal field calculation
%
% Calculation of focal field.
%
% See also Beam, BeamGauss.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Definition of beam

w0 = 5e-3;
Ex0 = 1;
Ey0 = 1i;
R = 10e-3;
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr)

P = 0.010;
bg = bg.normalize(P);

%% Calculation of focal field
 
dx = 10e-9;
dy = 10e-9;
[X,Y,Z] = meshgrid([-1000e-9:dx:1000e-9],[-1000e-9:dy:1000e-9],0e-6);

E = bg.focus(20e-3,X,Y,Z);
I = 0.5/PhysConst.Z0*sqrt(bg.er/bg.mr)*(E.*conj(E));

%% Plot of focal field

figure
title('Focal field')

surf(X*1e+9,Y*1e+9,Z*1e+9,I)

axis equal
shading interp
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

