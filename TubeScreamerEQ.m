[x, Fs_file] = audioread('guitarriff.wav');

if size(x,2) == 2
    x = mean(x,2);
end

Fs = 48000;
if Fs_file ~= Fs
    x = resample(x, Fs, Fs_file);
end

x = x / max(abs(x));

Drive = 0; % Value of 0-100
Tone = 0; % Value of 0-100
Volume = 100; % Value of 0-100

Drive  = min(max(Drive,  0), 100);
Tone   = min(max(Tone,   0), 100);
Volume = min(max(Volume, 0), 100);

%% Block 1
R  = 338000;
C  = 2e-8;
Av = 0.993;

% H1(s)
b1 = Av*R*C;
b0 = 0;
a1 = R*C;
a0 = 1;

num_1 = [b1 b0];
den_1 = [a1 a0];

[b_1,a_1] = bilinear(num_1, den_1, Fs);

x_hp = filter(b_1, a_1, x);
%{
[H,f] = freqz(b_1,a_1,4096,Fs);

figure;
subplot(2,1,1)
semilogx(f, 20*log10(abs(H)));
ylabel('Magnitude (dB)');
grid on;
title('Block 1 Frequency Response');

subplot(2,1,2)
semilogx(f, unwrap(angle(H)));
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
grid on;

figure;
impz(b_1, a_1);
title('Block 1 Impulse Response');

figure;
zplane(b_1, a_1);
title('Block 1 Pole-Zero Plot');
%}
%% Block 2
% High-pass Filter
R4 = 4.7e3;
C3 = 0.047e-6;

s = tf('s');

% H2(s)
H_hp = (s*R4*C3) / (1 + s*R4*C3);

[num, den] = tfdata(H_hp, 'v');
[b_dhp, a_dhp] = bilinear(num, den, Fs);

x_drivehp = filter(b_dhp, a_dhp, x_hp);

% Clipping Section
drive_alpha = Drive / 100; 

gain = 1 + 9*(drive_alpha^2); 
x_gain = gain .* x_drivehp; Vt = 0.3; 

x_nl = Vt * tanh(x_gain / Vt); 

x_clip = (1 - drive_alpha) * x_hp + drive_alpha * x_nl;

% Low-pass Filter
R6 = 51e3;
C4 = 51e-12;

% H3(s)
H_lp = 1 / (1 + s*R6*C4);

[num, den] = tfdata(H_lp, 'v');
[b_dlp, a_dlp] = bilinear(num, den, Fs);

x_drivetotal = filter(b_dlp, a_dlp, x_clip);

%% Block 3
T = Tone / 100;
Cz = 0.1e-6; 
Rz = 220;
Rf = 1000;
Rload = 10000;
Rpot = 20000;
Rl = max(T * Rpot, 1);
Rr = max((1 - T) * Rpot, 1);
wp1 = 1 / (Cz * (Rz + Rf));
wp2 = 1 / (Cz * Rload);
wz = 1 / (Cz * (Rz + Rl));
K = (Rload + Rf) / Rload;

% H4(s)
num = K * [1 wz];
den = conv([1 wp1], [1 wp2]);

[b_2, a_2] = bilinear(num, den, Fs);

peakGain = max(abs(H));
b_2 = b_2 / peakGain;

[H,f] = freqz(b_2,a_2,4096,Fs);

x_tone = filter(b_2, a_2, x_drivetotal);

rms_in  = sqrt(mean(x_drivetotal.^2));
rms_out = sqrt(mean(x_tone.^2));

if rms_out > 0
    x_tone = x_tone * (rms_in / rms_out);
end

volume_alpha = Volume / 100;
x_volume = volume_alpha * x_tone;
%{
figure;
subplot(2,1,1)
semilogx(f, 20*log10(abs(H)));
ylabel('Magnitude (dB)');
grid on;
title('Block 3 Frequency Response');

subplot(2,1,2)
semilogx(f, unwrap(angle(H)));
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
grid on;

figure;
impz(b_2, a_2);
title('Block 3 Impulse Response');

figure;
zplane(b_2, a_2);
title('Block 3 Pole-Zero Plot');
%}
%% Block 4
R  = 10000;
C  = 10e-6;

b1 = C*R;
b0 = 0;
a1 = C*R;
a0 = 1;

% H5(s)
num_3 = [b1 b0];
den_3 = [a1 a0];

[b_3,a_3] = bilinear(num_3, den_3, Fs);

x_output = filter(b_3, a_3, x_volume);

%{
[H,f] = freqz(b_3,a_3,4096,Fs);

figure;
subplot(2,1,1)
semilogx(f, 20*log10(abs(H)));
ylabel('Magnitude (dB)');
grid on;
title('Block 4 Frequency Response');

subplot(2,1,2)
semilogx(f, unwrap(angle(H)));
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
grid on;

figure;
impz(b_3, a_3);
title('Block 4 Impulse Response');

figure;
zplane(b_3, a_3);
title('Block 4 Pole-Zero Plot');
%}

%% Plotting Whole System Plots

b_total = conv(b_1, b_dhp);
b_total = conv(b_total, b_dlp);
b_total = conv(b_total, b_2);
b_total = conv(b_total, b_3);

a_total = conv(a_1, a_dhp);
a_total = conv(a_total, a_dlp);
a_total = conv(a_total, a_2);
a_total = conv(a_total, a_3);

[Htot, f] = freqz(b_total, a_total, 4096, Fs);

figure;
subplot(2,1,1)
semilogx(f, 20*log10(abs(Htot)));
ylabel('Magnitude (dB)');
grid on;
title('Total Linear System Frequency Response');

subplot(2,1,2)
semilogx(f, unwrap(angle(Htot)));
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
grid on;

figure;
impz(b_total, a_total);
title('Total Linear System Impulse Response');
xlim([0, 200]);

figure;
zplane(b_total, a_total);
title('Total Linear System Pole-Zero Plot');

%% Audio File
sound(x_output, Fs);

audiowrite('Drive100_Tone0_Volume80.wav', x_output, Fs);
