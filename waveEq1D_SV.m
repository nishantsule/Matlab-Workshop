%-------------------------------------
%%%% Solution to 1D wave equation
%-------------------------------------

%% Start by closing all open windows

close all  % close all existing figure windows

%% Display and user input 

% Enter code for displaying text and getting user input


%% if-else statements 

% Enter code to specify source parameters depending on the user input

imax = 800;  % total number of spatial grid points
% if the user enter l then run the following lines of code
%     % Light source
%     c = 3e8;  % speed of light in air [m/s]
%     cw = c / 1.33;  % speed of light in water [m/s]
%     f0 = 8e9;  % freq. of source [Hz]
%     amp = 0.8;  % source amplitude
%     ibd = round(imax / 3.5);  % location of material boundary
% else if the user enters s then run the following lines of code
%     % Sound source
%     c = 343;  % speed of sound in water [m/s]
%     cw = 1498;  % speed of sound in air [m/s]
%     f0 = 10000;  % freq. of source [Hz]
%     amp = 0.4;  % source amplitude
%     ibd = round(imax / 4);  % location of material boundary
% finally, if the user enters anything else display an error and stop
%     % Display error message 
%     error('You entered an invalid input.')


%% Initialize parameters

% Initiailize constant and array parameters for finite differences

lambda0 = min(c, cw) / f0;  % wavelength of source wave [m]
dx = lambda0 / 20;  % space grid step
dt = dx / max(c, cw);  % time grid step
nmax = round(0.5 * (imax * dx) / min(c, cw) / dt);  % number of time steps
w = 2 * pi * f0;  % angular frequency
tau = nmax * dt / 10;  % half width of source [s]
t0 = 3 * tau;  % time delay at source [s]
s1 = c * dt / dx;  % update coeff. for left half space
s2 = cw * dt / dx;  % update coeff. for right half space

% Define arrays
u = zeros(nmax, imax);  % initializing all fields to be zero
source = zeros(nmax, 1);  % initialize source array for all time steps

%% FDTD time-stepping for loop

% Enter code to update the variable 'u' in space and time

for n = 2 : nmax - 1
    
    % Define the source

    
    % Boundary conditions
    
    % At the left boundary: a source pulse propagating from left to right
%     u(n, 1) = source(n);
    % At the right boundary: all fields = 0 (reflecting boundary)
%     u(n, imax) = 0;
    
    % For loop to update u in space
    
    
    % Updating u in space without for loop (faster)
    
    
    % Plotting results
   
end