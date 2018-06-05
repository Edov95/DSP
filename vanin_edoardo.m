% clean up the environment
clear all; close all;

% read the .wav file
[signal_109,F] = audioread('signal_109.wav');

fprintf(' sampling frequency = %d Hz \n', F);
fprintf(' track duration = %5.2f s \n\n', length(signal_109)/F);

% the DFT
N = length(signal_109);     % analysis interval (relative to 10ms)
X = fft(signal_109(1:N));   % computation of the DFT of the whole signal
X_norm = X / N;             % normalization

% find the frequency of the carriers
locs = find(abs(X) >= (max(abs(X)) / 2));

A1 = abs(X_norm(locs(1))) + abs(X_norm(locs(4)));
A2 = abs(X_norm(locs(2))) + abs(X_norm(locs(3)));

freq = locs * F / N;

% debug simbol
%disp(freq);
%disp(A1);

%% Estranzione portante
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

%%%% Modulo e fase (?) dei filtri notch
[H_1, w] = freqz(b0*[1 -2 1], [1 -a1(1) -a2(1)], 'whole', 2048, F);
[H_2, w] = freqz(b0*[1 -2 1], [1 -a1(2) -a2(2)], 'whole', 2048, F);


%%%% Estraggo la portante
carry1 = filter(b0*[1 -2 1], [1 -a1(1) -a2(1)], signal_109);
carry2 = filter(b0*[1 -2 1], [1 -a1(2) -a2(2)], signal_109);

CARRY1 = fft(carry1(1:N)); %DFT of the carrier
CARRY1_NORM = CARRY1 / N;  %Normalization
CARRY2 = fft(carry2(1:N)); %DFT of the carrier
CARRY2_NORM = CARRY2 / N;  %Normalization


%sound(carry1, F);          %Suono le portanti
%sound(carry2, F);          %Suono le portanti

%%%%%%%%% Demodulation

left_demod  = signal_109 .* carry1 / A1;
right_demod = signal_109 .* carry2 / A2;

%sound(left_demod, F);
%% Clean & write

%%% highpass filter
f_3db_notch = 10;
delta_notch = pi*f_3db_notch/F;
r_notch = 1 - delta_notch;
a1 = -2*r_notch;
a2 = r_notch*r_notch;

b_highpass = [1 -2 1];
a_highpass = [1 a1 a2];
[H_highpass, w] = freqz(b_highpass, a_highpass, 'whole', 2048, F);

[b_lowpass, a_lowpass] = ellip(8, 5, 80, 2*4000 / F, 'low');
[H_lowpass, w] = freqz(b_lowpass, a_lowpass, 'whole', 2048, F);

left_lowpass = filter(10.*b_lowpass, a_lowpass, left_demod);
left = filter(10.*b_highpass, a_highpass, left_lowpass);

right_lowpass = filter(10.*b_lowpass, a_lowpass, right_demod);
right = filter(10.*b_highpass, a_highpass, right_lowpass);

y = zeros(N, 2);
y(:,1) = left;
y(:,2) = right;

audiowrite('signal_109_demod.wav', y, F);


%% plot figure 1
% Plot the magnitude
figure(1)
subplot(2,1,1);
plot(signal_109);
title('Signal in time');
xlabel('nT_c (s)'); ylabel('x(nT_c)')
subplot(2,1,2);                 % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
plot(f,20*log10(abs(X_norm)));
title('Magnitude (in dB) of the spectrum of the signal');
xlabel(' f (Hz)'); ylabel('|X_{norm}(f)|  (dB)');
axis([0 F -220 -20]);

%% plot figure 2
figure(2)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,N);              % frequency axis: 0---F Hz
subplot(2,1,1);
plot(f,20*log10(abs(CARRY1_NORM)));
title('Magnitude (in dB) of the spectrum of the carrier 1');
xlabel(' f (Hz)');
axis([0 F -220 -20]);
subplot(2,1,2);                 % Magnitude in dB (it is more meaningful) 
plot(f,20*log10(abs(CARRY2_NORM)));
title('Magnitude (in dB) of the spectrum of the carrier 2');
xlabel(' f (Hz)');
axis([0 F -220 -20]);

%% plot figure 4
figure(3)                       % Magnitude in dB (it is more meaningful)
f=linspace(0, F, 2048);              % frequency axis: 0---F Hz
subplot(2,1,1);
plot(f,20*log10(abs(H_highpass)));
title('Magnitude (in dB) of the spectrum of the high pass filter');
xlabel(' f (Hz)'); ylabel('|H_{highpass}(f)| (dB)');
maxy = max(20*log10(abs(H_highpass))) + 10; miny = maxy - 90;
axis([0 F miny maxy]);
subplot(2,1,2);
plot(f,angle(H_highpass));
title('Phase of the high pass filter');
xlabel(' f (Hz)'); ylabel('Phase(H_{highpass}(f)) (rad)');
xlim([0 F]);

%% plot figure 5
figure(4)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,2048);              % frequency axis: 0---F Hz
subplot(2,1,1);
plot(f,20*log10(abs(H_lowpass)));
title('Magnitude (in dB) of the spectrum of the low pass filter');
xlabel(' f (Hz)'); ylabel('|H_{lowpass}(f)|  (dB)');
maxy = max(20*log10(abs(H_lowpass))) + 10; miny = maxy - 90;
axis([0 F miny maxy]);
subplot(2,1,2);
plot(f,angle(H_lowpass));
title('Phase of the low pass filter');
xlabel(' f (Hz)'); ylabel('Phase(H_{lowpass}(f)) (rad)');
xlim([0 F]);

%% plot figure 6
figure(5)                       % Magnitude in dB (it is more meaningful)
f=linspace(0,F,2048);           % frequency axis: 0---F Hz
subplot(2,2,1);
plot(f,20*log10(abs(H_1)));
xlabel(' f (Hz)'); ylabel('|H_1(f)|  (dB)');
maxy = max(20*log10(abs(H_1))) + 10; miny = maxy - 90;
axis([0 F miny maxy]);

subplot(2,2,2);
plot(f,20*log10(abs(H_2)));
xlabel(' f (Hz)'); ylabel('|H_2(f)|  (dB)');
maxy = max(20*log10(abs(H_2))) + 10; miny = maxy - 90;
axis([0 F miny maxy]);

subplot(2,2,3);
plot(f,angle(H_1));
xlabel('f (Hz)'); ylabel('Phase(H_1(f)) (rad)');
xlim([0 F]);

subplot(2,2,4);
plot(f,angle(H_2));
xlabel(' f (Hz)'); ylabel('Phase(H_2(f)) (rad)');
xlim([0 F]);




