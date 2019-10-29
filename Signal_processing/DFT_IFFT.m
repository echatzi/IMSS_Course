%% Demo -- ifft
% =========================================================================
close all; clear; clc
m = 1;                          % system mass
omega = 5*pi;                   % natural freq.
zeta = 0.5;                     % damping ratio
omega_d = omega*sqrt(1-zeta^2); % natural freq. with damping
k =(omega^2)*m;                 % stiffness coefficient

dt = 0.005;                     % Sampling period 
t = 0:dt:20;                    % Time instants
F = cos(4.75*pi*t);             % Input force

% System impulse response -- time domain (convolution integral)
impulse = exp(-omega*zeta*t).*sin(omega_d*t);
output = dt*(1/m/(omega_d))*conv(impulse,F);

% System impulse response -- frequency domain (multiplying FRF by the FFT of the input signal)
output2 = dt*1/m/omega_d*ifft(fft(F).*fft(impulse));

% Figure -- comparison of the outputs calculated in time and freq. domain
figure(1)
subplot(211),plot(t,F(1:length(t)))
grid
title('Forcing Function')
subplot(212),plot(t,output(1:length(t)));hold on
plot(t,output2(1:length(t)),'--r')
legend({'Response calculated in time Domain','Response calculated in frequency domain'})
grid
title('System Response')
xlabel('Time (s)')



