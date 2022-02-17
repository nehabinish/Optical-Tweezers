% fieldpropagation.m - Field propagation along the z-axis
%
% Propagation of an electromagnetic field along the z-axis.
% Note: The units are not physical.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters
k = .5;

[x,y] = meshgrid([-10:.2:10],[-10:.2:10]);
z = pi;

E1.kx = 0;
E1.ky = .45;
E1.phi0 = -0.685;

E2.kx = -.03;
E2.ky = .30;
E2.phi0 = -1.253;

E3.kx = .02;
E3.ky = .15;
E3.phi0 = exp(i*E3.kx*x+i*E3.ky*y-i*1.497);

E4.kx = .35;
E4.ky = 0;
E4.phi0 = -1.122;

%% Initial fields

E1.field0 = exp(i*E1.kx*x+i*E1.ky*y+i*E1.phi0);
E2.field0 = exp(i*E2.kx*x+i*E2.ky*y+i*E2.phi0);
E3.field0 = exp(i*E3.kx*x+i*E3.ky*y+i*E3.phi0);
E4.field0 = exp(i*E4.kx*x+i*E4.ky*y+i*E4.phi0);

subplot(1,2,1)
hold on

Etot0 = E1.field0+E2.field0+E3.field0+E4.field0;
contourf(x,y,real(Etot0).^2,16)
shading flat
colormap([ones(100,1), [1:-0.01:0.01]', [1:-0.01:0.01]'])
axis equal tight

%% Propagated fields

E1.field1 = E1.field0*exp(i*sqrt(k^2-E1.kx^2-E1.ky^2)*z);
E2.field1 = E2.field0*exp(i*sqrt(k^2-E2.kx^2-E2.ky^2)*z);
E3.field1 = E3.field0*exp(i*sqrt(k^2-E3.kx^2-E3.ky^2)*z);
E4.field1 = E4.field0*exp(i*sqrt(k^2-E4.kx^2-E4.ky^2)*z);

subplot(1,2,2)
hold on

Etot1 = E1.field1+E2.field1+E3.field1+E4.field1;
contourf(x,y,real(Etot1).^2,16)
shading flat
colormap([ones(100,1), [1:-0.01:0.01]', [1:-0.01:0.01]'])
axis equal tight