%% Code based on the code offered by(c) dsprelated.com

%% Example 1: FFT of a 2 sine signal
clear all;close all;clc
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
figure(1);
subplot(2,1,1);
plot(n,x,'*k');
ni = [0:.1:N-1];     % Interpolated time axis
hold on;
plot(ni,A*cos(2*pi*ni*f1*fs+phi)+A*cos(2*pi*ni*f2*fs+phi),'-k'); grid off;
title('Sinusoid at 1/4 & 0.1 of the Sampling Rate');
xlabel('Time (samples)');
ylabel('Amplitude');
text(-8,1,'a)');
hold off;

% Plot spectral magnitude:
magX = abs(X);
fn = [0:1/N:1-1/N];  % Normalized frequency axis
subplot(2,1,2);
stem(fn,magX,'ok'); grid on;
xlabel('Normalized Frequency (cycles per sample)');
ylabel('Magnitude (Linear)');
text(-.11,40,'b)');
%% Correct the DFT using a Hanning Window
% Compute Hann window:
M = 31;         % Window length
nm = [0:M-1];   % time indices for window computation
% Hann window = "raised cosine", normalization (1/M)
% chosen to give spectral peak magnitude at 1/2:
w = (1/M) * (cos((pi/M)*(nm-(M-1)/2))).^2;

wzp = [w,zeros(1,N-M)]; % zero-pad out to the length of x
xw = x .* wzp;          % apply the window w to signal x

% Display real part of windowed signal and Hann window
figure
plot(n,wzp,'-k');grid on
title(['Hann Window']);
xlabel('Time (samples)'); ylabel('Amplitude');

% The Matlab for computing the DFT of the Hann-windowed complex sinusoid and plotting the results is listed below. 
% Compute the spectrum and its alternative forms:
Xw = fft(xw);              % FFT of windowed data
fn = [0:1.0/N:1-1.0/N];    % Normalized frequency axis

% Compute heavily interpolated versions for comparison:
Nzp = 16;                   % Zero-padding factor
Nfft = N*Nzp;               % Increased FFT size
xwi = [xw,zeros(1,Nfft-N)]; % New zero-padded FFT buffer
Xwi = fft(xwi);             % Compute interpolated spectrum
fni = [0:1.0/Nfft:1.0-1.0/Nfft]; % Normalized freq axis


figure;
plot(fn,abs(Xw),'*k'); hold on;
plot(fni,abs(Xwi),'-k'); hold off;
title('Spectral Magnitude');
xlabel('Normalized Frequency (cycles per sample))');
ylabel('Amplitude (linear)');


