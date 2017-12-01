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

%%%%%%% Calcolo parametri filtro passa banda per estrarre la portante

f_3_dB = 1;                 %Frequenza di taglio in Hz
theta_3_dB = 2*pi*f_3_dB/F; %Pulsazione di taglio

delta = theta_3_dB / 2;
r = 1 - delta;
b0 = delta;

a1 = zeros(2);
a2 = zeros(2);

for i = 1:2 
    a1(i) = 2*r*cos(2*pi*freq(i)/F);
    a2(i) = -r*r;
end
[H_1, w] = freqz(b0*[1 -2 1], [1 -a1(1) -a2(1)], 'whole', 2048, F);
[H_2, w] = freqz(b0*[1 -2 1], [1 -a1(2) -a2(2)], 'whole', 2048, F);

carry1 = filter(b0*[1 -2 1], [1 -a1(1) -a2(1)], signal_109);
carry2 = filter(b0*[1 -2 1], [1 -a1(2) -a2(2)], signal_109);

CARRY1 = fft(carry1(1:N));
CARRY1_NORM = CARRY1 / N;
CARRY2 = fft(carry2(1:N));
CARRY2_NORM = CARRY2 / N;


sound(carry1, F);
sound(carry2, F);





% Plot the magnitude
figure(1)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(X_norm)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_norm(f)|  (dB)');
axis([0 F -220 -20]);


figure(2)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(CARRY1_NORM)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_norm(f)|  (dB)');
axis([0 F -220 -20]);

figure(3)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(CARRY2_NORM)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_norm(f)|  (dB)');
axis([0 F -220 -20]);





