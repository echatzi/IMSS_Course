%% Code based on the code offered by(c) dsprelated.com

%% Example 1: FFT of a DFT-sinusoid
clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f = 0.25;            % Frequency (cycles/sample)
n = [0:N-1];         % Discrete time axis
x = A*cos(2*pi*n*f*fs+phi); % Sampled sinusoid
X = fft(x);          % Spectrum

% Plot time data:
figure(1);
subplot(2,3,1);
% subplot(3,1,1);

plot(n,x,'*k');
ni = [0:.1:N-1];     % Interpolated time axis
hold on;
plot(ni,A*cos(2*pi*ni*f*fs+phi),'-k'); grid off;
title('Sinusoid at 1/4 the Sampling Rate');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plot spectral magnitude:
figure(2)
subplot(2,3,1);
magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear)');
title('Sinusoid at 1/4 the Sampling Rate');


%% Example 2: FFT of a 2 sine signal
% clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f1 = 0.25;            % Frequency (cycles/sample) on a DFT line
f2 = 0.1;            % Frequency (cycles/sample) between two DFT lines

n = [0:N-1];         % Discrete time axis
x = A*cos(2*pi*n*f1*fs+phi)+A*cos(2*pi*n*f2*fs+phi); % Sampled 2-sine wave
X = fft(x);          % Spectrum

% Plot time data:
% figure(2);
figure(1);
subplot(2,3,2);

% subplot(2,1,1);
plot(n,x,'*k');
ni = [0:.1:N-1];     % Interpolated time axis
hold on;
plot(ni,A*cos(2*pi*ni*f1*fs+phi)+A*cos(2*pi*ni*f2*fs+phi),'-k'); grid off;
title('Sinusoid at 1/4 & 0.1');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plot spectral magnitude:
figure(2);
subplot(2,3,2);

magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear)');
title('Sinusoid at 1/4 & 0.1');

%% Example 3: FFT of a step function
% clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f1 = 0.25;            % Frequency (cycles/sample) on a DFT line
f2 = 0.1;            % Frequency (cycles/sample) between two DFT lines

%Construct a Step
for n=0:N-1
   if n>N/2-1
       x(n+1)=1;
   else
       x(n+1)=-1;
   end
end
n = [0:N-1];         % Discrete time axis

X = fft(x);          % Spectrum

% Plot time data:
figure(1);
subplot(2,3,3);
plot(n,x,'k');hold on;plot(n,x,'*k');
title('Step Function');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plot spectral magnitude:
figure(2);
subplot(2,3,3);
magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear)');
title('Step Function');

%% Example 4: FFT of an impulse function
% clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f1 = 0.25;            % Frequency (cycles/sample) on a DFT line
f2 = 0.1;            % Frequency (cycles/sample) between two DFT lines

%Construct an Impulse
for n=0:N-1
   if n==N/2-1
       x(n+1)=1;
   else
       x(n+1)=0;
   end
end
n = [0:N-1];         % Discrete time axis

X = fft(x);          % Spectrum

% Plot time data:
figure(1);
subplot(2,3,4);
plot(n,x,'k');hold on;plot(n,x,'*k');
title('Impulse Function');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plot spectral magnitude:
figure(2);
subplot(2,3,4);
magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
title('Impulse Function');
ylabel('Magnitude (Linear)');

%% Example 5: FFT of a Constant function
% clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f1 = 0.25;            % Frequency (cycles/sample) on a DFT line
f2 = 0.1;            % Frequency (cycles/sample) between two DFT lines

n = [0:N-1];         % Discrete time axis
x(n+1)=5;            % Construct a constant function

X = fft(x);          % Spectrum

% Plot time data:
figure(1);
subplot(2,3,5);
plot(n,x,'k');hold on;plot(n,x,'*k');
title('Constant Function');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plot spectral magnitude:
figure(2);
subplot(2,3,5);
magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (Linear)');
title('Constant Function');

%% Example 6: FFT of another DFT-sinusoid
% clear all;close all;clc
% Parameters:
N = 64;              % Must be a power of two
fs = 1;              % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f = 0.125;           % Frequency (cycles/sample)
n = [0:N-1];         % Discrete time axis
x = A*cos(2*pi*n*f*fs+phi); % Sampled sinusoid
X = fft(x);          % Spectrum

% Plot time data:
figure(1);
subplot(2,3,6);
plot(n,x,'*k');
ni = [0:.1:N-1];     % Interpolated time axis
hold on;
plot(ni,A*cos(2*pi*ni*f*fs+phi),'-k'); grid off;
% title('Sinusoid at .125Hz');
xlabel('Time (samples)');
ylabel('Amplitude');
hold off;

% Plofn(1)t spectral magnitude:
magX = abs(X);
fn = [0:1/N:1-1/N]*fs;  % frequency axis

figure(2);
subplot(2,3,6);
stem(fn,magX,'ok'); grid on;
xlabel('Frequency (Hz)');
title('Constant Function');
ylabel('Magnitude (Linear)');
xlim([0 .5])


