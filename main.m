% clc; close all; clear all;
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

%% Quantization of low power elements

%Fourier Transform
fft1 = fft(seg1);
fft2 = fft(seg2);
fft3 = fft(seg3);
N1 = length(seg1);

omega1khz = (-floor(N1/2):(N1-1-floor(N1/2)))*(fs/N1)/1000;

figure(1); plot(omega1khz, fftshift_norm_db(fft1));
title("FFT Magnitude (Segment 1)")
xlabel("Frequency (kHz)")
ylabel("|FFT|")

%% Remove >20 kHz
% sounds muffled
fft1Mid = floor(length(fft1)/2);
LPFcutoff = 20000 * N1/fs;
fft1LPF = fft1;
fft1LPF(LPFcutoff+2:fft1Mid) = 0;
fft1LPF(fft1Mid:end-LPFcutoff) = 0;
lpf1_out = ifft(fft1LPF);
% sound(seg1(1:10*fs), fs);
% pause();
% sound(lpf1_out(1:10*fs), fs);
figure(2); plot(omega1khz, fftshift_norm_db(fft1LPF)); title("lpf");

%% Quantization of low power elements
fft1Quant = fft1;
threshold = 0.00001 * max(abs(fft1));
lowIndex = abs(fft1) < threshold;
fft1Quant(lowIndex) = 0;

figure(3); plot(omega1khz, fftshift_norm_db(fft1Quant)); title("quantised");
seg1Out = ifft(fft1Quant);
error = (seg1Out - seg1);
figure(4)
plot((1:length(error(160:end-160)))/fs,error(160:end-160));
title( ...
    sprintf("Error in time domain of removing frequencies\nwith power under %.3f%% of max power", ...
    threshold*100/max(abs(fft1))) ...
    )
xlabel("Seconds")
ylabel("Absolute Error")
% sound(seg1_out(1:fs*5), fs);

%% Static thresholds
thresholds = [0.000001:0.000001:0.00005-0.000001,0.00005:0.00001:0.001];
[errors, coeffDiff, numCoeffDel] = test_static_thresholds(fft1, thresholds, seg1);
figure(5);
plot(thresholds*max(abs(fft1)), 20*log10(abs(errors)))
title("Normalised Error by Threshold Level")
xlabel("Lowest magnitude kept")
ylabel("Normalised Error (dB)")

figure(6);
plot(numCoeffDel, 20*log10(abs(errors)))
title("Normalised Error by Number of Coefficients Deleted")
xlabel("Number of Coefficients Deleted")
ylabel("Normalised Error (dB)")

figure(7);
plot(thresholds*100,coeffDiff);
title("Number of new coefficients deleted at current threshold compared to previous threshold")
xlabel("Thresholds (% of max magnitude of power)")
ylabel("No. of new coeff deleted")

%% test a threshold
fft1Quant = fft1;
threshold = 0.0001 * max(abs(fft1));
lowIndex62 = abs(fft1) < threshold;
fft1Quant(lowIndex62) = 0;
seg1Out = ifft(fft1Quant);
error62 = (seg1Out - seg1);
ssError = mean(abs(error62(160:end-160)))
sound(seg1Out(1:7*fs), fs);
% pause();
figure(8);
plot(omega1khz, fftshift_norm_db(fft1Quant));
title("threshold = 0.00062");
temp1 = fft1Quant;

%% Dynamic Thresholds and Smaller windows
% Base thresholds around peaks 
% Split signal into smaller overlapping blocks and apply the thresholds 
% for each block
winSize = 2048; % Size of block
overlap = 0.5; % Amount of overlap
maskedMod  = 0.5; % value to reduce masked frequencies (multiply)

seg1Out = window_signal(seg1, length(fft1), winSize, overlap, maskedMod);
sound(seg1Out(1:10*fs), fs);

figure(9);
plot(omega1khz, fftshift_norm_db(fft(seg1Out)));
title("dynamic threshold");