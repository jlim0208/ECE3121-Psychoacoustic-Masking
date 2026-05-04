% clc; close all; clear all;
[y, fs] = audioread('02Tchaikovsky_CapriccioItalienOp.45_cut.mp3');

%% Segment audio
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
N1= length(seg1);

omega1khz = (-floor(N1/2):(N1-1-floor(N1/2)))*(fs/N1)/1000;

function [output] = fftshiftnormdb(fft)
    fftshifted = fftshift(abs(fft));
    output = 20*log10(abs(fftshifted)/max(abs(fftshifted)));
end

figure(1); plot(omega1khz, fftshiftnormdb(fft1));
title("FFT Magnitude (Segment 1)")
xlabel("Frequency (kHz)")
ylabel("|FFT|")

%% Remove >20 kHz
fft1_half = floor(length(fft1)/2);
LPFcutoff = 80000;
fft1_LPF = fft1;
fft1_LPF(LPFcutoff+2:fft1_half) = 0;
fft1_LPF(fft1_half:end-LPFcutoff) = 0;
lpf1_out = ifft(fft1_LPF);
% sound(lpf1_out(1:5*fs), fs);

%Quantization of low power elements
fft1_quant = fft1;
threshold = 0.00001 * max(abs(fft1));
low_index = abs(fft1) < threshold;
fft1_quant(low_index) = 0;

figure(2); plot(omega1khz, fftshiftnormdb(fft1_quant)); title("quantised");
figure(4); plot(omega1khz, fftshiftnormdb(fft1_LPF)); title("lpf");
seg1_out = ifft(fft1_quant);
error = (seg1_out - seg1);
figure(3)
plot(error);
% sound(seg1_out(1:fs*5), fs);

%%
clear errors
thresholds = [0.000001:0.000001:0.00005-0.000001,0.00005:0.00001:0.001];
errors = [0];
difference = [0];
prev = zeros(1,length(fft1));
% errors = zeros(length(thresholds));
for i = thresholds
    fft1_quant = fft1;
    threshold = i * max(abs(fft1));
    low_index = abs(fft1) < threshold;
    fft1_quant(low_index) = 0;

    % ifft & error
    seg1_out = ifft(fft1_quant);
    error = (seg1_out - seg1);
    % Ignore 1st/last 160 values as there are edge artefacts
    av_error = mean(abs(error(160:end-160)));
    errors(end+1) = av_error;
    difference(end+1) = sum(low_index-prev);
    prev = low_index;
    % errors(floor(i/0.00005)+1) = av_error;
end
errors(1) = [];
difference(1) = [];

%%
figure(5);
plot(thresholds*100000,10*log10(abs(errors)))
figure(9);
plot(thresholds*100000,difference);

%% test the jump of error
fft1_quant = fft1;
threshold = 0.0001 * max(abs(fft1));
low_index62 = abs(fft1) < threshold;
fft1_quant(low_index62) = 0;
seg1_out = ifft(fft1_quant);
error62 = (seg1_out - seg1);
av_error = mean(abs(error62(160:end-160)))
sound(seg1_out(1:7*fs), fs);
% pause();
figure(6);
plot(omega1khz, fftshiftnormdb(fft1_quant));
title("threshold = 0.00062");
temp1 = fft1_quant;


fft1_quant = fft1;
threshold = 0.00063 * max(abs(fft1));
low_index63 = abs(fft1) < threshold;
fft1_quant(low_index63) = 0;
seg1_out = ifft(fft1_quant);
error63 = (seg1_out - seg1);
av_error = mean(abs(error63(160:end-160)))
% sound(seg1_out(1:7*fs), fs);

figure(7);
plot(omega1khz, fftshiftnormdb(fft1_quant));
title("threshold = 0.00063")
temp2 = fft1_quant;

%%
figure(8);
% temp1(isinf(temp1)) = 0;
% temp2(isinf(temp2)) = 0;
temp3 = temp1 - temp2;
temp3 = fftshiftnormdb(temp3);
% temp3(isinf(temp3))=0;

% plot(omega1khz, abs(temp3)/max(abs(temp3)));

% plot errors in subplot
subplot(2,1,1);
plot(error62(160:end-160));
title("error62");
subplot(2,1,2);
plot(error63(160:end-160));
title("error63");

%%
figure(10);
plot(seg1)