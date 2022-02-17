% focustm.m - Force on a sphere by a focused beam
%
% Calculation of the force produced by a focused beam on a spherical particle. 
% The sphere is kept at the origin of the coordinate system (i.e., its
% T-matrix is fixed), while the amplitudes of the focused beam are shifted.
% The results of a T-matrix calculation are compared to the exact
% analytical formula (Mie theory). 
%
% See also TMATRIX, TMATRIXSPHERE.

%   Author: Masoumeh Mousavi, Agnese Callegari, Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters and inizialization

% Beam
Ex0 = 1;  % x-component of electric field [V/m]
Ey0 = 0;  % y-component of electric field [V/m]
w0 = 5e-3;  % beam waist [m]
R = 5e-3;  % iris aperture [m]
lambda0 = 532e-9;  % vacuum wavelength [m]
k0 = 2*pi/lambda0;  % vacuum wavenumber [m^-1]
P = 0.010;  % power [W]

Nphi = 8;  % number of azimuathal divisions
Nr = 5;  % number of radial divisions

b = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);
b = b.normalize(P);

% focus
nm = 1.33;  % medium refractive index
NA = 1.30;  % objective numerical aperture
f = R*nm/NA;  % objective focal length [m]

fieldb = IncidentFieldFocusedBeam(b,f,nm);

% Mie particle
np = 1.50;  % refractive index
a = 100e-9;  % radius [m]

mie = MieParticle(nm,np,a,k0);
L = mie.lmax('formula','wiscombe');

tm = TMatrixSphere(mie,'Li',L,'Ls',L);

C = Point(-100e-9,-100e-9,-300e-9);  % Particle position [m]

%% Force calculation - T-matrix
Wi = fieldb.Wi(L,C);  % coefficients of the focused field
Ftm = tm.force(Wi)

%% Force calculation - Direct integration

% Parameters
er1 = 1.00;
er2 = nm^2;
mr1 = 1.00;
mr2 = 1.00;

n1 = sqrt(er1*mr1);
n2 = sqrt(er2*mr1);

% Integration parameters
radius = 4e-6;  % integration radius
N = 200;  % integration mesh
lambda = lambda0/nm;  % wavelength in the medium
QF = (2*pi*radius/N)/lambda;  % integration quality factor
disp(['Integration quality factory = ' num2str(QF) ' (should be < 1)'])

% Electromagnetic fields
[Theta,Phi,r] = meshgrid([pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N],radius);
[X,Y,Z] = Transform.Sph2Car(Theta,Phi,r);
X = X-C.X;
Y = Y-C.Y;
Z = Z-C.Z;

phi = b.phi;
dphi = phi(2,1)-phi(1,1);

theta = asin(b.r/f);
dr = b.r(1,2)-b.r(1,1);
R = dr*size(b.r,2);
dtheta = ones(size(b.r,1),1)*(asin([dr:dr:R]/f)-asin([0:dr:R-dr]/f));

k = 2*pi*sqrt(er2*mr2)/b.lambda0;
rho = sqrt(X.^2+Y.^2);
varphi = atan2(Y,X);

E = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
B = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
P = Point(X,Y,Z);
for n = 1:1:numel(theta)

    disp(['Calculation plane wave ' int2str(n) '/' int2str(numel(theta))])

    utheta = Point(cos(phi(n)).*cos(theta(n)), sin(phi(n)).*cos(theta(n)), -sin(theta(n)));
    utheta = utheta.tovector();
    uphi = Point(-sin(phi(n)), cos(phi(n)), zeros(size(phi(n))));
    uphi = uphi.tovector();

    % Radial component
    Pr = P.zrotation(-phi(n)).yrotation(-theta(n));
    [theta_r,phi_r,r_r] = Transform.Car2Sph(Pr.X,Pr.Y,Pr.Z);
    [Er_i,Br_i] = mie.incoming(theta_r,phi_r,r_r);
    [Er_s,Br_s] = mie.scattering(theta_r,phi_r,r_r,'L',L);
    Er_t = Er_i+Er_s;
    Er_t = Er_t.yrotation(theta(n)).zrotation(phi(n));
    Br_t = Br_i+Br_s;
    Br_t = Br_t.yrotation(theta(n)).zrotation(phi(n));

    % Azimuthal component    
    Pphi = P.zrotation(-phi(n)).yrotation(-theta(n)).zrotation(-pi/2);
    [theta_phi,phi_phi,r_phi] = Transform.Car2Sph(Pphi.X,Pphi.Y,Pphi.Z);
    [Ephi_i,Bphi_i] = mie.incoming(theta_phi,phi_phi,r_phi);
    [Ephi_s,Bphi_s] = mie.scattering(theta_phi,phi_phi,r_phi,'L',L);
    Ephi_t = Ephi_i+Ephi_s;
    Ephi_t = Ephi_t.zrotation(pi/2).yrotation(theta(n)).zrotation(phi(n));
    Bphi_t = Bphi_i+Bphi_s;
    Bphi_t = Bphi_t.zrotation(pi/2).yrotation(theta(n)).zrotation(phi(n));

    Etot = (b.Ephi(n)*Ephi_t + b.Er(n)*Er_t)*(sqrt(n1/n2)*sqrt(cos(theta(n))));
    Btot = (b.Ephi(n)*Bphi_t + b.Er(n)*Br_t)*(sqrt(n1/n2)*sqrt(cos(theta(n))));

    shift = exp(1i * (k*sin(theta(n)).*cos(phi(n))*C.X+k*sin(theta(n)).*sin(phi(n))*C.Y+k*cos(theta(n))*C.Z) );
    E = E + ComplexVector(E.X,E.Y,E.Z, ...
        Etot.Vx*shift*sin(theta(n))*dphi*dtheta(n), ...
        Etot.Vy*shift*sin(theta(n))*dphi*dtheta(n), ...
        Etot.Vz*shift*sin(theta(n))*dphi*dtheta(n) ...
        );
    B = B + ComplexVector(B.X,B.Y,B.Z, ...
        Btot.Vx*shift*sin(theta(n))*dphi*dtheta(n), ...
        Btot.Vy*shift*sin(theta(n))*dphi*dtheta(n), ...
        Btot.Vz*shift*sin(theta(n))*dphi*dtheta(n) ...
        );
end
E = 1i*k*f*exp(-1i*k*f)/(2*pi)*E;
B = 1i*k*f*exp(-1i*k*f)/(2*pi)*B;

% Forces
[ux,uy,uz] = Transform.Sph2CarVector(Theta,Phi,zeros(size(r)),zeros(size(r)),ones(size(r)));
ur = Vector(E.X,E.Y,E.Z,ux,uy,uz);

dOmega = sin(Theta) * radius * (Theta(1,2)-Theta(1,1)) * radius * (Phi(2,1)-Phi(1,1));

ExEr = ComplexVector(E.X,E.Y,E.Z, ...
    E.Vx .* (conj(E).*ur), ...
    E.Vy .* (conj(E).*ur), ...
    E.Vz .* (conj(E).*ur) ...
    );
BxBr = ComplexVector(B.X,B.Y,B.Z, ...
    B.Vx .* (conj(B).*ur), ...
    B.Vy .* (conj(B).*ur), ...
    B.Vz .* (conj(B).*ur) ...
    );
df = .5*PhysConst.e0*mie.nm^2 .* ( ...
    ( ExEr + (PhysConst.c0./mie.nm)^2*BxBr ) .* dOmega ...
    - .5*( E.norm().^2 + (PhysConst.c0./mie.nm)^2*B.norm().^2 ) .* ur .* dOmega ...
    ) ;

Fint = Vector(0,0,0,real( sum(sum(df.Vx)) ),real( sum(sum(df.Vy)) ),real( sum(sum(df.Vz)) ) )

%% Final report
Ftm
Fint