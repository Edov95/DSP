% clean up the environment
clear all; close all;

% read the .wav file
[signal_109,F] = audioread('signal_109.wav');

fprintf(' sampling frequency = %d Hz \n', F);
fprintf(' track duration = %5.2f s \n\n', length(signal_109)/F);

% the DFT
N = length(signal_109);     % analysis interval (relative to 10ms)
XX = fft(signal_109(1:N));  % computation of the DFT of the whole signal
XX = XX / N;                % normalization

% find the frequency of the carriers
freq = find(abs(XX) >= (max(abs(XX)) / 2));

A1 = (abs(XX(freq(1))) + abs(XX(freq(4)))) / 2;
A2 = (abs(XX(freq(2))) + abs(XX(freq(3)))) / 2;

freq = freq * F / N;

% debug simbol
disp(freq);

% Plot the magnitude
figure(1)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(XX)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|Y(f)|  (dB)');
maxy = max(20*log10(abs(XX))); miny=maxy-90;
axis([0 F miny maxy]);


Fstop1 = (freq(1) - 5);
Fpass1 = (freq(1) + 5);

b = firpm(500, [0 Fstop1-10 Fstop1 Fpass1 Fpass1+10 F]/F,[0 0 1 1 0 0]);

[H, W] = freqz(b, 1, length(XX));

disp(length(H));
disp(length(XX));

%Estraggo la portante
CARR1 = XX .* H;

carr1 = fft(CARR1(1:N)); %portante sinusoidale

%sound(abs(carr1), F);

left = signal_109 .* carr1 / A1; %demodulo

b4000 = firpm(500, [0 5 10 4000 4500 F]/F,[0 0 1 1 0 0]);

left_filtered = filter(b4000, 1, left);

sound(abs(left_filtered), F);










