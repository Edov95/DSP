% clean up the environment
clear all; close all;

% read the .wav file
[signal_109,F] = audioread('signal_109.wav');

fprintf(' sampling frequency = %d Hz \n', F);
fprintf(' track duration = %5.2f s \n\n', length(signal_109)/F);

% Listen the original signal
%fprintf('You are listening to an audio signal \n');
%sound(signal_109,F);

% the DFT
N = 1*F;                      % analysis interval (relative to 10ms)
XX = fft(signal_109(6500:6500 + N - 1));  % computation of the DFT of 10 ms (6500 is the starting point of the analysis window)T

% Computation of the notch filter coefficients
b1 = -2*cos(2*pi*17600/F); 
b = [1 b1 1];
delta = pi*5000/F; r=1-delta;
a2=r^2; a1 = -2*r*cos(2*pi*17600/F);
a= [1 a1 a2];

y=filter(b, a, signal_109);             % filter the signal

Y = fft(y(6500:6500 + N - 1));

% Plot the magnitude
figure(1)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz

subplot(2,1,1);
plot(f,20*log10(abs(XX)));
title('Magnitude (in dB) of the spectrum of the portion of the signal');
xlabel(' f (Hz)'); ylabel('|Y(f)|  (dB)');
maxy = max(20*log10(abs(XX))); miny=maxy-90;
axis([0 F miny  maxy]);
subplot(2,1,2);
plot(f,20*log10(abs(Y)));
title('Magnitude (in dB) of the spectrum of the portion of the signal');
xlabel(' f (Hz)'); ylabel('|Y(f)|  (dB)');
maxy = max(20*log10(abs(Y))); miny=maxy-90;
axis([0 F miny  maxy]);


sound(y,F);




