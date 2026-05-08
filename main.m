% close all; clear all; clc;
[y, fs] = audioread('02Tchaikovsky_CapriccioItalienOp.45_cut.mp3');

%% Segment audio
% segment 1: single brass
% segment 2: multiple brass
% segment 3: strings & some brass
cut12 = 18;
cut23 = 50;
seg1 = y(1:cut12*fs);
seg2 = y(cut12*fs+1:cut23*fs);
seg3 = y(cut23*fs+1:end);
% sound(seg3,fs)

%% Compute FFTs

%Fourier Transform
fft1 = fft(seg1);
fft2 = fft(seg2);
fft3 = fft(seg3);
N1 = length(seg1);

omega1khz = (-floor(N1/2):(N1-1-floor(N1/2)))*(fs/N1)/1000;

figure(1); plot(omega1khz, fftshift_norm_db(fft1));
title("FFT of Signal")
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

%% Remove >20 kHz
fft1Mid = floor(length(fft1)/2);
LPFcutoff = 20000 * N1/fs;
fft1LPF = fft1;
fft1LPF(LPFcutoff+2:fft1Mid) = 0;
fft1LPF(fft1Mid:end-LPFcutoff) = 0;
lpf1_out = ifft(fft1LPF);
% sound(seg1(1:10*fs), fs);
% pause();
% sound(lpf1_out(1:10*fs), fs);

figure(2); 
plot(omega1khz, fftshift_norm_db(fft1LPF)); 
title("Frequency domain of Low Pass Filter (Wc=20kHz) of signal");
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

%% Quantisation of low power elements
fftQuant1 = fft1;
threshold = 0.00001 * max(abs(fft1));
lowIndex = abs(fft1) < threshold;
fftQuant1(lowIndex) = 0;

figure(3); plot(omega1khz, fftshift_norm_db(fftQuant1)); 
title(sprintf("Filtered Frequency Domain \n(low-power components under 0.001%% of max removed)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

figure(4)
seg1Out1 = ifft(fftQuant1);
error1 = (seg1Out1 - seg1);
plot((1:length(error1(160:end-160)))/fs,error1(160:end-160));
title(sprintf("Error in Time Domain \n(low-power components under 0.001%% of max removed)"));
xlabel("Seconds")
ylabel("Absolute Error")
% sound(seg1_out(1:fs*5), fs);

fftQuant5 = fft1;
threshold = 0.0005 * max(abs(fft1));
lowIndex5 = abs(fft1) < threshold;
fftQuant5(lowIndex5) = 0;
seg1Out5 = ifft(fftQuant5);
error5 = (seg1Out5 - seg1);
ssError5 = sum(error5(160:end-160) .^ 2);
fprintf("Static sum squared error: %f\n", ssError5);
% sound(seg1Out(1:7*fs), fs);
% pause();
figure(5);
plot(omega1khz, fftshift_norm_db(fftQuant5));
title(sprintf("Filtered Frequency Domain \n(low-power components under 0.05%% of max removed)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

%% Static thresholds
thresholds = [0.000001:0.000001:0.00005-0.000001,0.00005:0.00001:0.001];
[errors, coeffDiff, numCoeffDel] = test_thresholds(fft1, thresholds * max(abs(fft1)), seg1, @static_masking);
figure(6);
plot(thresholds*100, 20*log10(abs(errors)))
title("Normalised Error by Threshold Level")
xlabel("Thresholds (% of max magnitude of power)")
ylabel("Normalised Error (dB)")

figure(7);
plot(numCoeffDel, 20*log10(abs(errors)))
title("Normalised Error by Number of Coefficients Deleted")
xlabel("Number of Coefficients Deleted")
ylabel("Normalised Error (dB)")

figure(8);
plot(thresholds*100,coeffDiff);
title("Number of new coefficients deleted at current threshold compared to previous threshold")
xlabel("Thresholds (% of max magnitude of power)")
ylabel("No. of new coeff deleted")

%% Dynamic Thresholds and Smaller windows
% Base thresholds around peaks 
% Split signal into smaller overlapping blocks and apply the thresholds 
% for each block
winSize = 2048; % Size of block
overlap = 0.5; % Amount of overlap
maskedMod  = 0.5; % value to reduce masked frequencies (multiply)

[seg1OutDyn, ~] = window_signal(seg1, fft1, winSize, overlap, maskedMod, @dynamic_masking);
sound(seg1OutDyn(1:10*fs), fs);

figure(9);
plot(omega1khz, fftshift_norm_db(fft(seg1OutDyn)));
title(sprintf("Filtered Frequency Domain \n(low-power components near peaks are halved)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

% Plot SSError
errorDyn = (seg1OutDyn - seg1);
ssErrorDyn = mean(errorDyn(160:end-160) .^ 2);
fprintf("Dynamic sum squared error: %f\n", ssErrorDyn);

%% Plot a snippet of the filtered signal in time domain
figure(10);
subplot(2,1,1);
tIdx = 5*fs:5.1*fs; % Time vector for plotting
t = tIdx/fs;
plot(t, seg1(tIdx));
title('Original Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, seg1OutDyn(tIdx));
title('Filtered Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

%% Test multiple masked modifiers for dynamic thresholds
maskedMods = 0:0.1:1;
[errorsMM, ~] = test_thresholds(fft1, maskedMods, seg1, @(fftSig, maskedMult) window_signal(seg1, fftSig, winSize, overlap, maskedMult, @dynamic_masking));

figure(11);
plot(maskedMods, 20*log10(abs(errorsMM)))
title("Normalised Error by Masking Modifier")
xlabel("Modifier Value")
ylabel("Normalised Error (dB)")

%% Combine plots for comparison
% Combine frequency plots for comparison
figure(12);
plot(omega1khz, fftshift_norm_db(fft1LPF), 'DisplayName', 'LPF (20kHz)');
hold on;
plot(omega1khz, fftshift_norm_db(fftQuant1), 'DisplayName', 'Quantised powers below 0.001%');
plot(omega1khz, fftshift_norm_db(fftQuant5), 'DisplayName', 'Quantised powers below 0.05%');
plot(omega1khz, fftshift_norm_db(fft(seg1OutDyn)), 'DisplayName', 'Quantised (halved) powers near peaks');
legend('Location','south');
title('Comparison of Frequency Domains with Different Filtering');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
hold off;