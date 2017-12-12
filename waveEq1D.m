%% Solution to 1D wave equation

disp('Do you want to simulate a light or a sound pulse?'); 
sim_flag = input('Enter l for light and s for sound): ', 's');

if sim_flag == 'l'
    % Light source
    c = 3e8;  % speed of light in air [m/s]
    cw = c / 2;  % speed of light in water [m/s]
    f0 = 6e9;  % freq. of source [Hz]
    amp = 1;
elseif sim_flag == 's'
    % Sound source
    c = 1498;  % speed of sound in water [m/s]
    cw = 343;  % speed of sound in air [m/s]
    f0 = 20000;  % freq. of source [Hz]
    amp = 1;
else
    % Display error message 
    error('You entered an invalid input.')
end

% FD parameters
imax = 500;  % total number of spatial grid points
dx = lambda0 / 20;  % space grid step
dt = dx / max(c, cw);  % time grid step
nmax = round(0.5 * (imax * dx) / min(c, cw) / dt);
w = 2 * pi * f0;  % angular frequency
lambda0 = min(c, cw) / f0;  % wavelength of source wave [m]
tau = nmax * dt / 50;  % half width of source [s]
t0 = 3 * tau;  % time delay at source [s]
s1 = c * dt / dx;  % update coeff. for left half space
s2 = cw * dt / dx;  % update coeff. for right half space

% Define array u
u = zeros(nmax, imax);  % initializing all fields to be zero

% FDTD time stepping loop
for n = 2 : nmax - 1
    
    % Boundary conditions
    % At the left boundary: a source pulse propagating from left to right
    u(n, 1) = amp * sin(w * (dt * n - t0)) ...
              .* exp(-((dt * n - t0).^2) / tau^2);
    % At the right boundary: all fields = 0 (reflecting boundary)
    u(n, imax) = 0;
    
    % Space update for loop
    for i = 2 : imax - 1
        
        if i < imax / 2
            u(n + 1, i) = s1^2 * (u(n, i + 1) ...
                          - (2 * u(n, i)) + u(n, i - 1)) ...
                          + (2 * u(n, i)) - u(n - 1, i);
        else
            u(n + 1, i) = s2^2 * (u(n, i + 1) ...
                          - (2 * u(n, i)) + u(n, i - 1)) ...
                          + (2 * u(n, i)) - u(n - 1, i);
        end
        
    end
    % Space update without for loop (faster)
%     u(n+1, 2 : imax / 2) = s1^2 * (u(n, 3 : imax / 2 + 1) - ...
%                          (2 * u(n, 2 : imax / 2)) + ...
%                          u(n, 1 : imax / 2 - 1))...
%                          + (2 * u(n, 2 : imax / 2)) ...
%                          - u(n-1, 2 : imax / 2);
%     u(n+1, imax / 2 + 1 : imax - 1) = s2^2 * (u(n, imax / 2 + 2 : imax) - ...
%                          (2 * u(n, imax / 2 + 1 : imax - 1)) + ...
%                          u(n, imax / 2 : imax - 2))...
%                          + (2 * u(n, imax / 2 + 1 : imax - 1)) ...
%                          - u(n-1, imax / 2 + 1 : imax - 1);
    
    % Plotting results
    plot(u(n + 1, 1:imax))
    hold on
    plot([imax / 2, imax / 2], [-1, 1], '--k')
    hold off
    axis([1, imax, -1, 1])
    xlabel('x'), ylabel('u')
    T = ['Time step = ', num2str(n + 1), ' of ', num2str(nmax)];
    title(T);
    pause(0.01)
end