% bessel.m - Simulation of a Bessel beam
%
% Simulation of a Bessel beam.
% Note: the units are not physical.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters
phi = 2*pi*[.001:.001:1];
A = ones(size(phi));

L = 64;

kp = 1;

[x,y] = meshgrid([-20:.1:20],[-20:.1:20]);

%% Calculation
kx = kp*cos(phi);
ky = kp*sin(phi);

E = 0*x;
for i = 1:1:length(phi)
    E = E+A(i)*exp(1i*kx(i)*x+1i*ky(i)*y);
end
I = E.*conj(E);
a = angle(E);
a(a<-3.14) = 3.14;

%% Plots

figure
set(gcf,'color','w')
axes('position',[0 0 1 1])
contourf(x,y,I,L)
colormap(ones(L,3)-[0:1/(L-1):1]'*[0 1 1])
shading flat
axis([-20 20 -20 20])
axis off

figure
set(gcf,'color','w')
axes('position',[0 0 1 1])
contourf(x,y,a,L)
shading flat
colormap(ones(L,3)-[0:1/(L-1):1]'*[.5 .5 .5]-.25)
axis([-20 20 -20 20])
axis off