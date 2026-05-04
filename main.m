% clc; close all; clear all;
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

function [output] = fftshiftnormdb(fft)
fftshifted = fftshift(abs(fft));
output = 20*log10(abs(fftshifted)/max(abs(fftshifted)));
end

figure(1); plot(omega1, fftshiftnormdb(fft1));
title("FFT Magnitude (Segment 1)")
xlabel("Frequency (Hz)")
ylabel("|FFT|")

%% Remove >20 kHz
fft1_half = floor(length(fft1)/2);
LPFcutoff = 80000;
fft1_LPF = fft1;
fft1_LPF(LPFcutoff+2:fft1_half) = 0;
fft1_LPF(fft1_half:end-LPFcutoff) = 0;
% sound(ifft(fft1_LPF), fs);

%Quantization of low power elements
fft1_quant = fft1;
threshold = 0.00001 * max(abs(fft1));
low_index = abs(fft1) < threshold;
fft1_quant(low_index) = 0;

figure(2); plot(omega1, fftshiftnormdb(fft1_quant)); title("quantised");
figure(4); plot(omega1, fftshiftnormdb(fft1_LPF)); title("lpf");
seg1_out = ifft(fft1_quant);
error = (seg1_out - seg1);
figure(3)
plot(error);
sound(seg1_out(1:fs*5), fs);


% for i 