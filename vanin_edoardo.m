% clean up the environment
clear all; close all;

% read the .wav file
[signal_109,F] = audioread('signal_109.wav');

fprintf(' sampling frequency = %d Hz \n', F);
fprintf(' track duration = %5.2f s \n\n', length(signal_109)/F);

% the DFT
N = length(signal_109);     % analysis interval (relative to 10ms)
X = fft(signal_109(1:N));  % computation of the DFT of the whole signal
X_norm = X / N;                % normalization

% find the frequency of the carriers
locs = find(abs(X) >= (max(abs(X)) / 2));

A1 = abs(X_norm(locs(1))) + abs(X_norm(locs(4)));
A2 = abs(X_norm(locs(2))) + abs(X_norm(locs(3)));

freq = locs * F / N;

% debug simbol
disp(freq);
disp(A1);
disp(A2);

% Plot the magnitude
figure(1)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(X_norm)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_norm(f)|  (dB)');
axis([0 F -220 -20]);


%%%%%%% Calcolo parametri filtro passa banda per estrarre la portante

f_3_dB = 1;                 %Frequenza di taglio in Hz
theta_3_dB = 2*pi*f_3_dB/F; %Pulsazione di taglio

delta = theta_3_dB / 2;
r = 1 - delta;
b0 = delta;
a1 = 2*r*cos(2*pi*freq(1)/F);
a2 = -r*r;

[H_1, w] = freqz(b0*[1 -2 1], [1 -a1 -a2], 'whole', 2048, F);

carry1 = filter(b0*[1 -2 1], [1 -a1 -a2], signal_109);

CARRY1 = fft(carry1(1:N));

%sound(carry1, F);


figure(2)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(CARRY1)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_norm(f)|  (dB)');
axis([0 F -220 -20]);







