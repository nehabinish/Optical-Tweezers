% dipoleforces.m - Dipole forces near a focused optical beam
%
% Dipole forces near a focusef optical beam.
%
% See also InducedDipole, BeamGauss, EFieldFocus.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters
dx = 5e-9;
dy = 5e-9;
dz = 5e-9;
Lx = 305e-9;
Ly = 305e-9;
Lz = 605e-9;

L = 32;

%% Dipole

ep = 2.25;
a = 10e-9;
lambda0 = 632e-9;
alpharc = InducedDipole.polarizability('radiative correction',a,ep,'lambda0',lambda0);
id = InducedDipole(alpharc,lambda0);

%% Beam

w0 = 5e-3;
Ex0 = 1;
Ey0 = 0;
R = 5e-3;
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',632e-9);
bg = bg.normalize(0.010);

%% Focused field

nm = 1.33;
NA = 1.20;
f = nm*R/NA;

ef = EFieldFocus(bg,f)

%% Field
subplot(1,2,1)

[X,Y,Z] = meshgrid([-Lx:dx:Lx],[-Ly:dy:Ly],0);
E = bg.focus(f,X,Y,Z,'er',1.33^2);
I = .5*PhysConst.c0/nm*(nm^2*PhysConst.e0)*norm(E).^2;

contourf(X*1e+6,Y*1e+6,I,L)
colormap(ones(L,3)-[0:1/(L-1):1]'*[0 1 1])
shading flat
axis equal off
axis([-.40 .40 -.40 .40])

%% Forces
subplot(1,2,2)

[X,Y,Z] = meshgrid([-Lx:10*dx:Lx],[-Ly:10*dy:Ly],0);
r = Point(X,Y,Z);
[F,Fgrad,Fscat,Fsc] = id.force(r,ef);

h = maxis2d([-.37 .39],[-.37 .37], ...
    'color', [.7 .7 .7], ...
    'xlim',[-.40 .40],'xlabel','$x$', 'xticks', [], ...
    'ylim',[-.40 .40],'ylabel','$y$', 'yticks', []);
hold on
F = .06*F./max(max(abs(norm(F))));
arrow2d(X*1e+6, Y*1e+6, real(X*1e+6+F.Vx), real(Y*1e+6+F.Vy), ...
    'stemwidth',.001,'headlength',.01,'headnode',.01,'headwidth',.005,'linewidth',.05);
axis equal