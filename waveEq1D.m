%% Solution to 1D wave equation

imax = 1000;
nmax = 500;

% Light source
c = 3e8;  % speed of light in free space [m/s]
v = c / 2;  % speed of light in lossless material [m/s]
f0 = 3e9;  % freq. of source [Hz]
w = 2 * pi * f0;  % angular frequency
lambda0 = c / f0;  % wavelength of source wave [m]
tau = 0.4e-9;  % half width of source [s]
t0 = 1.5e-9;  % time delay at source [s]

% FD parameters
dx = lambda0 / 20;  % space grid step
dt = dx * S / c;  % time grid step
s1 = c * dt / dx;  % update coeff. for left half space
s2 = v * dt / dx;  % update coeff. for right half space

% Define array u
u = zeros(nmax, imax);  % initializing all fields to be zero

% FDTD time stepping loop
for n = 2 : nmax - 1
    
    % Boundary conditions
    % At the left boundary: a source pulse propagating from left to right
    u(n,1) = sin(w*(dt*n-t0)).*exp(-((dt*n-t0).^2)/tau^2);
    % At the right boundary: all fields = 0 (reflecting boundary)
    u(n,imax) = 0;
    
    % Space update for loop
%     for i = 2 : imax - 1
%         
%         if i < imax / 2
%             u(n+1,i) = s1^2 * (u(n, i+1) - (2 * u(n, i)) + u(n, i-1))...
%                         + (2 * u(n, i)) - u(n-1, i);
%         else
%             u(n+1, i) = s2^2 * (u(n, i+1) - (2 * u(n, i)) + u(n, i-1))...
%                          + (2 * u(n, i)) - u(n-1, i);
%         end
%         
%     end
    % Space update without for loop (faster)
    u(n+1, 2 : imax / 2) = s1^2 * (u(n, 3 : imax / 2 + 1) - ...
                         (2 * u(n, 2 : imax / 2)) + ...
                         u(n, 1 : imax / 2 - 1))...
                         + (2 * u(n, 2 : imax / 2)) ...
                         - u(n-1, 2 : imax / 2);
    u(n+1, imax / 2 + 1 : imax - 1) = s1^2 * (u(n, imax / 2 + 2 : imax) - ...
                         (2 * u(n, imax / 2 + 1 : imax - 1)) + ...
                         u(n, imax / 2 : imax - 2))...
                         + (2 * u(n, imax / 2 + 1 : imax - 1)) ...
                         - u(n-1, imax / 2 + 1 : imax - 1);
    
    % Plotting results
    plot(u(n, 1:imax))
    axis([0, 500, -1, 1])
    xlabel('x'), ylabel('u')
    T = ['N_{\lambda}=', num2str(Nl),'  S=', num2str(S), ...
        '  Time steps=', num2str(nmax)];
    title(T);
    pause(0.01)
    
end