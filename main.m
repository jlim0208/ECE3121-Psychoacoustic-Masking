clc; close all; clear all;
[y, fs] = audioread('02Tchaikovsky_CapriccioItalienOp.45_cut.mp3');

%% Segment audio
cut12 = 18;
cut23 = 50;
seg1 = y(1:cut12*fs);
seg2 = y(cut12*fs+1:cut23*fs);
seg3 = y(cut23*fs+1:end);
sound(seg3,fs)

%% Quantization of low power elements

%Fourier Transform
fft1 = fft(seg1);
fft2 = fft(seg2);
fft3 = fft(seg3);
N1= length(seg1);

omega1 = (-floor(N1/2):(N1-1-floor(N1/2)))*(fs/N1);
figure(1); plot(omega1, fftshift(abs(fft1)));
title("FFT Magnitude (Segment 1)")
xlabel("Frequency (Hz)")
ylabel("|FFT|")

%Quantization of low power elements
fft1_quant = fft1;
threshold = 0.05 * max(abs(fft1));
low_index = abs(fft1) < threshold;
fft1_quant(low_index) = round(fft1(low_index) / 100) * 100;

figure(2); plot(omega1, fftshift(abs(fft1_quant)));
seg1_out = ifft(fft1_quant);
sound(seg1_out, fs);