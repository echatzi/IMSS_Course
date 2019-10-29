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
xlabel('Normalized Frequency (cycles per sample)');
ylabel('Amplitude (linear)');

%% Compare different Window options on a single sine signal
clear all;close all;clc
% Parameters:
N = 64*2;              % Must be a power of two
fs = 1;               % Set sampling rate to 1 Hz
A = 1;               % Sinusoidal amplitude
phi = 0;             % Sinusoidal phase
f2 = 0.1;            % Frequency (cycles/sample) between two DFT lines

n = [0:N-1];         % Discrete time axis
x = A*cos(2*pi*n*f2*fs+phi); % Sampled 2-sine wave
X = fft(x);          % Spectrum

% Compute Hann window using standard MATLAB functions:
M = N/2-1;         % Window length
wH=hann(M);
wzpH = [wH',zeros(1,N-M)]; % zero-pad out to the length of x
xwH = x .* wzpH;          % apply the window w to signal x

% Display real part of windowed signal and Hann window
figure (100)
subplot(2,2,1)
plot(n,wzpH,'-k');grid on
title(['Hann Window']);
xlabel('Time (samples)'); ylabel('Amplitude');

% The Matlab for computing the DFT of the Hann-windowed complex sinusoid and plotting the results is listed below. 
% Compute the spectrum and its alternative forms:
XwH = fft(xwH);              % FFT of windowed data
fn = [0:1.0/N:1-1.0/N];    % Normalized frequency axis

% Compute heavily interpolated versions for comparison:
Nzp = 16;                   % Zero-padding factor
Nfft = N*Nzp;               % Increased FFT size
xwiH = [xwH,zeros(1,Nfft-N)]; % New zero-padded FFT buffer
XwiH = fft(xwiH);             % Compute interpolated spectrum
fni = [0:1.0/Nfft:1.0-1.0/Nfft]; % Normalized freq axis

figure(200);
subplot(2,2,1)
plot(fn,abs(XwH),'*k'); hold on;
plot(fni,abs(XwiH),'-k'); hold off;
title('Hanning Window');
xlabel('Normalized Frequency (cycles per sample)');
ylabel('Amplitude (linear)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Rectangular window using standard MATLAB functions:
wR=rectwin(M);

wzpR = [wR',zeros(1,N-M)]; % zero-pad out to the length of x
xwR = x .* wzpR;          % apply the window w to signal x

% Display real part of windowed signal and Hann window
figure (100)
subplot(2,2,2)
plot(n,wzpR,'-k');grid on
title(['Rectangular Window']);
xlabel('Time (samples)'); ylabel('Amplitude');

% The Matlab for computing the DFT of the Hann-windowed complex sinusoid and plotting the results is listed below. 
% Compute the spectrum and its alternative forms:
XwR = fft(xwR);              % FFT of windowed data
fn = [0:1.0/N:1-1.0/N];    % Normalized frequency axis

% Compute heavily interpolated versions for comparison:
Nzp = 16;                   % Zero-padding factor
Nfft = N*Nzp;               % Increased FFT size
xwiR = [xwR,zeros(1,Nfft-N)]; % New zero-padded FFT buffer
XwiR = fft(xwiR);             % Compute interpolated spectrum
fni = [0:1.0/Nfft:1.0-1.0/Nfft]; % Normalized freq axis

figure(200);
subplot(2,2,2)
plot(fn,abs(XwR),'*k'); hold on;
plot(fni,abs(XwiR),'-k'); hold off;
title('Rectangular Window');
xlabel('Normalized Frequency (cycles per sample))');
ylabel('Amplitude (linear)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Flattop window using standard MATLAB functions:
wF=flattopwin(M);

wzpF = [wF',zeros(1,N-M)]; % zero-pad out to the length of x
xwF = x .* wzpF;          % apply the window w to signal x

% Display real part of windowed signal and Hann window
figure (100)
subplot(2,2,3)
plot(n,wzpF,'-k');grid on
title(['Flattop Window']);
xlabel('Time (samples)'); ylabel('Amplitude');

% The Matlab for computing the DFT of the Hann-windowed complex sinusoid and plotting the results is listed below. 
% Compute the spectrum and its alternative forms:
XwF = fft(xwF);              % FFT of windowed data
fn = [0:1.0/N:1-1.0/N];    % Normalized frequency axis

% Compute heavily interpolated versions for comparison:
Nzp = 16;                   % Zero-padding factor
Nfft = N*Nzp;               % Increased FFT size
xwiF = [xwF,zeros(1,Nfft-N)]; % New zero-padded FFT buffer
XwiF = fft(xwiF);             % Compute interpolated spectrum
fni = [0:1.0/Nfft:1.0-1.0/Nfft]; % Normalized freq axis

figure(200);
subplot(2,2,3)
plot(fn,abs(XwF),'*k'); hold on;
plot(fni,abs(XwiF),'-k'); hold off;
title('Flattop Window');
xlabel('Normalized Frequency (cycles per sample)');
ylabel('Amplitude (linear)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Blackman window using standard MATLAB functions:
wB=blackman(M);

wzpB = [wB',zeros(1,N-M)]; % zero-pad out to the length of x
xwB = x .* wzpB;          % apply the window w to signal x

% Display real part of windowed signal and Hann window
figure (100)
subplot(2,2,4)
plot(n,wzpB,'-k');grid on
title(['Blackman Window']);
xlabel('Time (samples)'); ylabel('Amplitude');

% The Matlab for computing the DFT of the Hann-windowed complex sinusoid and plotting the results is listed below. 
% Compute the spectrum and its alternative forms:
XwB = fft(xwB);              % FFT of windowed data
fn = [0:1.0/N:1-1.0/N];    % Normalized frequency axis

% Compute heavily interpolated versions for comparison:
Nzp = 16;                   % Zero-padding factor
Nfft = N*Nzp;               % Increased FFT size
xwiB = [xwB,zeros(1,Nfft-N)]; % New zero-padded FFT buffer
XwiB = fft(xwiB);             % Compute interpolated spectrum
fni = [0:1.0/Nfft:1.0-1.0/Nfft]; % Normalized freq axis

figure(200);
subplot(2,2,4)
plot(fn,abs(XwB),'*k'); hold on;
plot(fni,abs(XwiB),'-k'); hold off;
title('Blackman Window');
xlabel('Normalized Frequency (cycles per sample)');
ylabel('Amplitude (linear)');







