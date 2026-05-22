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
N2 = length(seg2);
N3 = length(seg3);

omega1khz = (-floor(N1/2):(N1-1-floor(N1/2)))*(fs/N1)/1000;
omega2khz = (-floor(N2/2):(N2-1-floor(N2/2)))*(fs/N2)/1000;
omega3khz = (-floor(N3/2):(N3-1-floor(N3/2)))*(fs/N3)/1000;

figure(1); plot(omega1khz, fftshift_norm_db(fft1));
title("FFT of Signal")
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

%% Remove >20 kHz
fft1Mid = floor(length(fft1)/2);
LPFcutoff = 10000 * N1/fs;
fft1LPF = fft1;
fft1LPF(LPFcutoff+1:fft1Mid) = 0;
fft1LPF(fft1Mid:end-LPFcutoff+1) = 0;
lpf1Out = ifft(fft1LPF);
% sound(seg1(1:10*fs), fs);
% pause();
% sound(lpf1Out(1:10*fs), fs);

figure(2); 
plot(omega1khz, fftshift_norm_db(fft1LPF)); 
title("Frequency domain of Low Pass Filter (Wc=10kHz) of signal");
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

% Find SSError
errorLPF = (lpf1Out - seg1);
ssErrorLPF = mean(errorLPF(160:end-160) .^ 2);
fprintf("LPF sum squared error (sig1): %e\n", ssErrorLPF);

%% Quantisation of low power elements
[fftQuant1, ~] = static_masking(fft1, 0.00001*max(abs(fft1)));

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

[fftQuant5, s1o5DelIndexes] = static_masking(fft1, 0.0005*max(abs(fft1)));
seg1Out5 = ifft(fftQuant5);
error5 = (seg1Out5 - seg1);
ssError5 = sum(error5(160:end-160) .^ 2);
fprintf("Static sum squared error (sig1): %e\n", ssError5);
% sound(seg1Out(1:7*fs), fs);
% pause();
figure(5);
plot(omega1khz, fftshift_norm_db(fftQuant5));
title(sprintf("Filtered Frequency Domain \n(low-power components under 0.05%% of max removed)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

%% Static thresholds
thresholds = [0.000001:0.000001:0.00005-0.000001,0.00005:0.00001:0.001];
[errorsSig1, coeffDiffSig1, numCoeffDelSig1] = test_thresholds(fft1, thresholds * max(abs(fft1)), seg1, @static_masking);
[errorsSig2, coeffDiffSig2, numCoeffDelSig2] = test_thresholds(fft2, thresholds * max(abs(fft2)), seg2, @static_masking);
[errorsSig3, coeffDiffSig3, numCoeffDelSig3] = test_thresholds(fft3, thresholds * max(abs(fft3)), seg3, @static_masking);

%% Plot error against thresholds
figure(6);
plot(thresholds*100, 20*log10(abs(errorsSig1)),'DisplayName', 'Single Brass')
hold on;
plot(thresholds*100, 20*log10(abs(errorsSig2)),'DisplayName', 'Multiple Brass')
plot(thresholds*100, 20*log10(abs(errorsSig3)),'DisplayName', 'Strings and Some Brass')
hold off;
title("Normalised Error by Threshold Level")
xlabel("Thresholds (% of max magnitude of power)")
ylabel("Normalised Error (dB)")
legend('Location','southeast')

figure(7);
plot(numCoeffDelSig1, 20*log10(abs(errorsSig1)),'DisplayName', 'Single Brass')
hold on;
plot(numCoeffDelSig1, 20*log10(abs(errorsSig2)),'DisplayName', 'Multiple Brass')
plot(numCoeffDelSig1, 20*log10(abs(errorsSig3)),'DisplayName', 'Strings and Some Brass')
hold off;
title("Normalised Error by Number of Coefficients Deleted")
xlabel("Number of Coefficients Deleted")
ylabel("Normalised Error (dB)")
legend('Location','southeast')

figure(8);
plot(thresholds*100, log10(coeffDiffSig1),'DisplayName', 'Single Brass');
hold on;
plot(thresholds*100, log10(coeffDiffSig2),'DisplayName', 'Multiple Brass');
plot(thresholds*100, log10(coeffDiffSig3),'DisplayName', 'Strings and Some Brass');
hold off;
title("Difference in Coefficients Deleted between Thresholds")
xlabel("Thresholds (% of max magnitude of power)")
ylabel("No. of new coeff deleted (log10)")
legend('Location','southeast')

%% Dynamic Thresholds and Smaller windows
% Base thresholds around peaks 
% Split signal into smaller overlapping blocks and apply the thresholds 
% for each block
winSize = 2048; % Size of block
overlap = 0.5; % Amount of overlap
maskGain  = 0.5; % value to reduce masked frequencies (multiply)

[fft1OutDyn, s1DynDelIndexes] = window_signal(seg1, fft1, winSize, overlap, maskGain, @dynamic_masking);
seg1OutDyn = ifft(fft1OutDyn);
sound(seg1OutDyn(1:10*fs), fs);

figure(9);
plot(omega1khz, fftshift_norm_db(fft(seg1OutDyn)));
title(sprintf("Filtered Frequency Domain \n(low-power components near peaks are halved)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

% Find SSError
errorDyn = (seg1OutDyn - seg1);
ssErrorDyn = mean(errorDyn(160:end-160) .^ 2);
fprintf("Dynamic sum squared error (sig1): %e\n", ssErrorDyn);

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

%% Test multiple attenuation factors for dynamic thresholds
maskGains = 0:0.1:1;
[errorsMM1, ~] = test_thresholds(fft1, maskGains, seg1, @(fftSig, maskedMult) window_signal(seg1, fftSig, winSize, overlap, maskedMult, @dynamic_masking));
[errorsMM2, ~] = test_thresholds(fft2, maskGains, seg2, @(fftSig, maskedMult) window_signal(seg2, fftSig, winSize, overlap, maskedMult, @dynamic_masking));
[errorsMM3, ~] = test_thresholds(fft3, maskGains, seg3, @(fftSig, maskedMult) window_signal(seg3, fftSig, winSize, overlap, maskedMult, @dynamic_masking));

%% Plot error against attenuation factors
close(findobj('Type', 'figure', 'Number', 11));
figure(11);
yyaxis right
plot(1-maskGains, abs(errorsMM1),'DisplayName', 'Single Brass')
yyaxis left
plot(1-maskGains, abs(errorsMM2),'DisplayName', 'Multiple Brass')
hold on
plot(1-maskGains, abs(errorsMM3),'DisplayName', 'Strings and Some Brass')
hold off
title("Error by Masking Attenuation Factor")
xlabel("Attenuation Factor (dB) of Masked Bins")
ylabel("Sum of Squared Error")
legend('Location','northwest')

%% Test multiple decays for dynamic thresholds
% Test how error changes against an increasing masking range
decaysList = num2cell(2:2:40);
for i = 1:length(decaysList)
    decaysList(i) = {[4 1] ./ decaysList{i}};
end
[errorsMD1, ~] = test_thresholds(fft1, decaysList, seg1, @(fftSig, decays) window_signal(seg1, fftSig, winSize, overlap, maskGain, @(x,y) dynamic_masking(x,y,decays)));
[errorsMD2, ~] = test_thresholds(fft2, decaysList, seg2, @(fftSig, decays) window_signal(seg2, fftSig, winSize, overlap, maskGain, @(x,y) dynamic_masking(x,y,decays)));
[errorsMD3, ~] = test_thresholds(fft3, decaysList, seg3, @(fftSig, decays) window_signal(seg3, fftSig, winSize, overlap, maskGain, @(x,y) dynamic_masking(x,y,decays)));

%% Plot error against decays
close(findobj('Type', 'figure', 'Number', 12));
figure(12);
yyaxis right
plot(2:2:40, abs(errorsMD1),'DisplayName', 'Single Brass')
yyaxis left
plot(2:2:40, abs(errorsMD2),'DisplayName', 'Multiple Brass')
hold on
plot(2:2:40, abs(errorsMD3),'DisplayName', 'Strings and Some Brass')
title("Error against Masking Range of Different Instruments")
xlabel(sprintf("Forward Masking Range Factor \n(Factor is Multiplied to the Range (bins) where Forward Decay is 1 dB/bin)"))
ylabel("Sum of Squared Error")
legend('Location','southeast')
hold off

%% Combine plots for comparison
% Combine frequency plots for comparison
figure(13);
hexCols = rgb2hex(orderedcolors("gem"));
ogPlt = plot(omega1khz, fftshift_norm_db(fft1), 'Color', 'k', 'DisplayName', 'Unfiltered');
hold on;
qntPlt1 = plot(omega1khz, fftshift_norm_db(fftQuant1), 'Color', hexCols(2), 'DisplayName', 'Quantised powers below 0.001%');
lpfPlt  = plot(omega1khz, fftshift_norm_db(fft1LPF), 'Color', hexCols(1), 'DisplayName', 'LPF (10kHz)');
qntPlt5 = plot(omega1khz, fftshift_norm_db(fftQuant5), 'Color', hexCols(3), 'DisplayName', 'Quantised powers below 0.05%');
dynPlt  = plot(omega1khz, fftshift_norm_db(fft(seg1OutDyn)), 'Color', hexCols(4), 'DisplayName', 'Quantised (halved) powers near peaks');
legend([ogPlt, lpfPlt, qntPlt1, qntPlt5, dynPlt], 'Location','south');
title('Comparison of Frequency Domains with Different Filtering');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
hold off;

%% Test dynamic for multiple brass
winSize = 2048; % Size of block
overlap = 0.5; % Amount of overlap
maskGain  = 0.5; % value to reduce masked frequencies (multiply)

[fft2OutDyn, ~] = window_signal(seg2, fft2, winSize, overlap, maskGain, @dynamic_masking);
seg2OutDyn = ifft(fft2OutDyn);
sound(seg2OutDyn(1:10*fs), fs);

figure(14);
plot(omega2khz, fftshift_norm_db(fft(seg2OutDyn)));
title(sprintf("Filtered Frequency Domain of multiple brass instruments \n(low-power components near peaks are halved)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

% Find SSError
errorDyn2 = (seg2OutDyn - seg2);
ssErrorDyn2 = mean(errorDyn2(160:end-160) .^ 2);
fprintf("Dynamic sum squared error (sig2): %e\n", ssErrorDyn2);

%% Test dynamic for strings and some brass
[fft3OutDyn, ~] = window_signal(seg3, fft3, winSize, overlap, maskGain, @dynamic_masking);
seg3OutDyn = ifft(fft3OutDyn);
sound(seg3OutDyn(1:10*fs), fs);

figure(15);
plot(omega3khz, fftshift_norm_db(fft(seg3OutDyn)));
title(sprintf("Filtered Frequency Domain of string and some brass instruments\n(low-power components near peaks are halved)"));
xlabel("Frequency (kHz)")
ylabel("Magnitude (dB)")

% Find SSError
errorDyn3 = (seg3OutDyn - seg3);
ssErrorDyn3 = mean(errorDyn3(160:end-160) .^ 2);
fprintf("Dynamic sum squared error (sig3): %e\n", ssErrorDyn3);

%% Save Figures/Outputs (generated by gemini, modified)
% Get list of all open figure handles
figList = findobj(allchild(0), 'flat', 'Type', 'figure');
mkdir('results')

% Loop through and save each figure
for iFig = 1:length(figList)
    figHandle = figList(iFig);

    figNum = num2str(get(figHandle, 'Number'));

    % Save the figure
    exportgraphics(figHandle, ['results/Figure ' figNum '.png']); 
end

fileID = fopen("results\output.txt", 'w');
fprintf(fileID, "LPF sum squared error (sig1): %e\n", ssErrorLPF);
fprintf(fileID, "Static sum squared error (sig1): %e\n", ssError5);
fprintf(fileID, "Dynamic sum squared error (sig1): %e\n", ssErrorDyn);
fprintf(fileID, "Dynamic sum squared error (sig2): %e\n", ssErrorDyn2);
fprintf(fileID, "Dynamic sum squared error (sig3): %e", ssErrorDyn3);
fclose(fileID);

%% Write audio
audiowrite("Seg1_normal.mp3",seg1,fs);
audiowrite("Seg1_lpf.mp3",lpf1Out,fs);
audiowrite("Seg1_static_mask.mp3",seg1Out5,fs);
audiowrite("Seg1_dynamic_mask.mp3",seg1OutDyn,fs);

%% Calculate compression ratio
fprintf("Compression ratio of static masking: %f\n", sum(s1o5DelIndexes)/N1);
fprintf("Compression ratio of dynamic masking: %f\n", mean(s1DynDelIndexes/winSize));