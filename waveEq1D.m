%%%% Solution to 1D wave equation

close all  % close all existing figure windows

disp('Do you want to simulate a light or a sound pulse?'); 
sim_flag = input('Enter l for light or s for sound): ', 's');

imax = 800;  % total number of spatial grid points
if sim_flag == 'l'
    % Light source
    c = 3e8;  % speed of light in air [m/s]
    cw = c / 1.33;  % speed of light in water [m/s]
    f0 = 8e9;  % freq. of source [Hz]
    amp = 0.8;  % source amplitude
    ibd = round(imax / 3.5);  % location of material boundary
elseif sim_flag == 's'
    % Sound source
    c = 343;  % speed of sound in water [m/s]
    cw = 1498;  % speed of sound in air [m/s]
    f0 = 10000;  % freq. of source [Hz]
    amp = 0.4;  % source amplitude
    ibd = round(imax / 4);  % location of material boundary
else
    % Display error message 
    error('You entered an invalid input.')
end

% FD parameters
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

% FDTD time stepping loop
for n = 2 : nmax - 1
    
    % Define the source
    source(n) = amp * sin(w * (dt * n - t0)) ...
              .* exp(-((dt * n - t0).^2) / tau^2);
    
    % Boundary conditions
    % At the left boundary: a source pulse propagating from left to right
    u(n, 1) = source(n);
    % At the right boundary: all fields = 0 (reflecting boundary)
    u(n, imax) = 0;
    
    % For loop to update u in space
    for i = 2 : imax - 1
        
        if i < ibd
            u(n + 1, i) = s1^2 * (u(n, i + 1) ...
                          - (2 * u(n, i)) + u(n, i - 1)) ...
                          + (2 * u(n, i)) - u(n - 1, i);
        else
            u(n + 1, i) = s2^2 * (u(n, i + 1) ...
                          - (2 * u(n, i)) + u(n, i - 1)) ...
                          + (2 * u(n, i)) - u(n - 1, i);
        end
        
    end
    % Updating u in space without for loop (faster)
%     u(n+1, 2 : ibd) = s1^2 * (u(n, 3 : ibd + 1) - ...
%                          (2 * u(n, 2 : ibd)) + ...
%                          u(n, 1 : ibd - 1))...
%                          + (2 * u(n, 2 : ibd)) ...
%                          - u(n-1, 2 : ibd);
%     u(n+1, ibd + 1 : imax - 1) = s2^2 * (u(n, ibd + 2 : imax) - ...
%                          (2 * u(n, ibd + 1 : imax - 1)) + ...
%                          u(n, ibd : imax - 2))...
%                          + (2 * u(n, ibd + 1 : imax - 1)) ...
%                          - u(n-1, ibd + 1 : imax - 1);
    
    % Plotting results
    plot(u(n + 1, 1:imax))
    hold on
    plot([ibd, ibd], [-1, 1], '--k');
    hold off
    axis([1, imax, -1, 1]);
    xlabel('x'), ylabel('u');
    T = ['Air                        ','Time step = ', num2str(n + 1), ...
        ' of ', num2str(nmax), '                                  Water'];
    title(T);
    pause(0.01);
end