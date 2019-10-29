%% Demo 3 (pseudorandom signal created by sines)
% =========================================================================
clear; close all;clc;
N = 2000;                   % Number of samples
x = zeros(N,1);             % Memory allocation for vector x
Ts = 0.1;                   % Sampling period
Fs = 1/Ts;                  % Sampling frequency
t = 0:Ts:(N-1)*Ts;          % Time sequence
F = Fs*(0:N-1)/(N-1);       % Frequency grid
omega = 0.0:0.1:5;          % Sinus function frequencies
A = randn(length(omega),1); % Sinus function amplitudes

% p sinus waves of different amplitudes and frequencies added one by one to
% vector x
for p = 1:length(omega)     
    x = x + A(p)*sin(omega(p)*t*2*pi)';
    % Plot the time history and the corresponding FFT of x 
    if p<=5 || p == length(omega)
        figure(1)
        subplot(211),plot(t,x)
        xlabel('Time (s)')
        ylabel('x(t)')
        subplot(212),plot(F,abs(fft(x,N)))
        xlabel('Frequency (Hz)')
        ylabel('|F(\omega)|')
        pause
    end
end



%% Demo 3 -- ifft
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


%% Demo 3 -- Aliasing I
% =========================================================================
close all
N = 20000;                      % Number of samples
t = 0:1e-3:1e-3*(N-1);          % Fine time grid (for the approximation of the continuous curve) 
x = sin(2*pi*t);                % 1 Hz sine        

% Plot the sinus wave 
figure(1)
plot(t,x)   
axis tight
xlabel('Time (s)')
ylabel('Sinusoidal signal')
title('Sampling problems')
hold on 
Tsample = 1:200:N;                          % Sampling time vector
Ts = t(201);                                % Sampling period
Fs = 1/Ts;                                  % Sampling frequency
disp(['Sampling time: Ts = ',num2str(Ts),' s,  ','Sampling frequency fs = ',num2str(Fs),' (Hz)'])
pause(1)
plot(t(Tsample),x(Tsample),'--or')          % Plot sample points 
xlim([0 5])
pause(1)
plot(t,sin(10*pi*t),'-.g')                  % Plot sine wave of 6 Hz 
legend({'Original sine wave (1 Hz)','sample points (5Hz)','sine wave (5 Hz)'})
pause(3)
% Note: The original signal (blue line) may be considered as continuous 
%       However, once sampled there is no way to distinguish it from
%       signals with much higher frequency

% FFT of the sampled signal
N = length(Tsample);                        % Number of sampled points 
NFFT = 2^nextpow2(N);                       % Calculating the min power p with 2^p > N
Y = fft(x(Tsample),NFFT)/N;                 % FFT calculation
f = Fs/2*linspace(0,1,NFFT/2+1);            % Frequency points for the calculated FFT 

% Plot single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
grid on
% Note: In this case the frequency content of the sampled signal is below 
%       the Nyquist frequency (fs/2) and thus is correctly shown in the
%       amplitude spectrum at 1 Hz

%% Demo 3 -- Aliasing II
% =========================================================================
close all
N = 150000;                     % Number of samples
t = 0:1e-3:1e-3*(N-1);          % Fine time grid (for the approximation of the continuous curve) 
x = sin(2*pi*t);                % 1 Hz sine      

% Plot the sinus wave
figure(1)
plot(t,x)
axis tight
xlabel('Time (s)')
ylabel('Sinusoidal signal')
title('Sampling problems')
hold on 

Tsample = 1:801:N;                          % Sampling time vector
Ts = t(801);                                % Sampling period
Fs = 1/Ts;                                  % Sampling frequency
disp(['Sampling time: Ts = ',num2str(Ts),' s,  ','Sampling frequency fs = ',num2str(Fs),' (Hz)'])

plot(t(Tsample),x(Tsample),'--or')          % Plot sample points 
xlim([0 5])

% FFT of the sampled signal
N = length(Tsample);                        % Number of sampled points 
NFFT = 2^nextpow2(N);                       % Calculating the min power p with 2^p > N
Y = fft(x(Tsample),NFFT)/N;                 % FFT calculation
f = Fs/2*linspace(0,1,NFFT/2+1);            % Frequency points for the calculated FFT 

% Plot single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
grid on
axis tight
% Note: In this case the continuous signal is sampled with 1.25 Hz and thus  
%       the original freq. is shown as shown aliased at 0.25 Hz.

%% Demo 3 -- Filtering demo
% Read and play a wav audio file
[funky, f] = wavread('funky.wav');
wavplay(funky, f);

% Plot a portion of the waveform 
figure(200)
subplot(2,1,1), plot(funky), title('Entire waveform');
smallRange = 100000:100200;
subplot(2,1,2), plot(smallRange, funky(smallRange)), title('100 milliseconds');

% Plot the spectrogram -that is the frequency content over time
figure(400)
specgram(funky, 512, f);
% and of a portion of it
close all;
subplot(4,1,1), plot(funky), title('Entire waveform');
subplot(4,1,2),specgram(funky, 512, f);
subplot(4,1,3), plot(funky(100000:110000)), axis('tight');
title('Part of Audio')
subplot(4,1,4), specgram(funky(100000:110000),128,f);

%% Demo 3 -- Low-pass filter

% We design a Butterworth 10th order low-pass filter to supress frequencies
% higher than 600Hz (Sampling freq. = 22050 Hz). 
fNorm = 600 / (f/2);                % Cutoff freq. given normalized wrt the Nyquist freq.
[numfl,denfl] = butter(10, fNorm, 'low');     % Calculate filter's TF 
funkyLow = filtfilt(numf, denf, funky);     % Filter out the signal
wavplay(funkyLow, f);                       % Listen the result


% Filter's FRF
[MagL,PhaseL,WL] = dbode(numfl,denfl,1/f);

%% Demo 3 -- High-pass filter
% We design a Butterworth 10th order high-pass filter to supress frequencies
% lower than 5kHz (Sampling freq. = 22050 Hz). 
fNorm = 5000 / (f/2);                % Cutoff freq. given normalized wrt the Nyquist freq.
[numfh,denfh] = butter(10, fNorm, 'high');     % Calculate filter's TF 
funkyHigh = filtfilt(numfh, denfh, funky);     % Filter out the signal
wavplay(funkyHigh, f);                       % Listen the result

% Filter's FRF
[MagH,PhaseH,WH] = dbode(numfh,denfh,1/f);

%% Demo 3 -- Pass-Band filter
% We design a Butterworth 10th order band-pass filter to supress frequencies
% lower than 2kHz and higher than 2500 (Sampling freq. = 22050 Hz). 
fNorm = [2000/(f/2), 2500/(f/2)];            % Cutoff freq. given normalized wrt the Nyquist freq.
[numfb,denfb] = butter(10, fNorm, 'bandpass');   % Calculate filter's TF 
funkyBand = filtfilt(numfb, denfb, funky);   % Filter out the signal
wavplay(funkyBand, f);                       % Listen the result

% Filter's FRF
[MagB,PhaseB,WB] = dbode(numfb,denfb,1/f);

%% Demo 3 -- Plot Filter's FRF
figure(600)
subplot(211),plot(WL/2/pi,20*log(MagL),'-k')
hold on
plot(WH/2/pi,20*log(MagH),'--b')
plot(WB/2/pi,20*log(MagB),'-.r')
xlabel('Frequency (Hz)')
ylabel('Magnitude (db)')
grid on
legend('Lowpass','Highpass','Bandpass')
subplot(212),plot(WL/2/pi,PhaseL,'-k')
hold on
plot(WH/2/pi,PhaseH,'--b')
plot(WB/2/pi,PhaseB,'-.r')
xlabel('Frequency (Hz)')
ylabel('Phase')
grid on