%% Solution to 1D wave equation

imax = 1000;
nmax = 500;

% Light source
c = 3e8;  % speed of light in free space [m/s]
v = c/2;  % speed of light in lossless material [m/s]
f0 = 3e9;  % freq. of source [Hz]
w = 2*pi*f0;  % angular frequency
lambda0 = c/f0;  % wavelength of source wave [m]
tau = 0.4e-9;  % half width of source [s]
t0 = 1.5e-9;  % time delay at source [s]

% FD parameters
dx = lambda0/Nl;    %space grid step
dt = dx*S/c;        %time grid step
s1 = c*dt/dx;       %Courant number for left half space
s2 = v*dt/dx;       %Courant number for right half space