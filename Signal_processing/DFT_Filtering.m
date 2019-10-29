clear all;close all;clc
%% Filtering demo
% Read and play a wav audio file
[funky, f] = audioread('beatit.mp3');
funky = funky(3.1e6:3.3e6);
playerObj=audioplayer(funky, f);
play(playerObj);

% Plot a portion of the waveform 
figure(200)
subplot(2,1,1), plot(funky), title('Entire waveform');
smallRange = 1:200;
subplot(2,1,2), plot(smallRange, funky(smallRange)), title('100 milliseconds');

% Plot the spectrogram -that is the frequency content over time
figure(400)
specgram(funky, 512, f);
% and of a portion of it
close all;
subplot(4,1,1), plot(funky), title('Entire waveform');
subplot(4,1,2),specgram(funky, 512, f);
subplot(4,1,3), plot(funky(1:10000)), axis('tight');
title('Part of Audio')
subplot(4,1,4), specgram(funky(1:10000),512,f);
ylim([0 1000])

%% Demo 3 -- Low-pass filter

% We design a Butterworth 10th order low-pass filter to supress frequencies
% higher than 600Hz. 
fNorm = 100 / (f/2);                % Cutoff freq. given normalized wrt the Nyquist freq.
[numfl,denfl] = butter(4, fNorm, 'low');     % Calculate filter's TF 
funkyLow = filtfilt(numfl, denfl, funky);     % Filter out the signal
playerObj=audioplayer(10000*funkyLow, f);
play(playerObj);

% Filter's FRF
[MagL,PhaseL,WL] = dbode(numfl,denfl,1/f);

%% High-pass filter
% We design a Butterworth 10th order high-pass filter to supress frequencies
% lower than 5kHz (Sampling freq. = 22050 Hz). 
fNorm = 4000 / (f/2);                % Cutoff freq. given normalized wrt the Nyquist freq.
[numfh,denfh] = butter(10, fNorm, 'high');     % Calculate filter's TF 
funkyHigh = filtfilt(numfh, denfh, funky);     % Filter out the signal
playerObj=audioplayer(10*funkyHigh, f);
play(playerObj);

% Filter's FRF
[MagH,PhaseH,WH] = dbode(numfh,denfh,1/f);

%% Pass-Band filter
% We design a Butterworth 10th order band-pass filter to supress frequencies
fNorm = [500/(f/2),900/(f/2)];            % Cutoff freq. given normalized wrt the Nyquist freq.
[numfb,denfb] = butter(4, fNorm, 'bandpass');   % Calculate filter's TF 
funkyBand = filtfilt(numfb, denfb, funky);   % Filter out the signal
playerObj=audioplayer(10*funkyBand, f);
play(playerObj);

% Filter's FRF
[MagB,PhaseB,WB] = dbode(numfb,denfb,1/f);

%% Plot Filter's FRF
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