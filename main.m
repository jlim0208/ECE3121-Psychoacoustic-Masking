% clc; close all; clear all;
[y, fs] = audioread('02Tchaikovsky_CapriccioItalienOp.45_cut.mp3');

%% Shift fft to 0-center, normalise and use dB scale
function [output] = fftshiftNormdB(fft)
    fftShifted = fftshift(abs(fft));
    output = 20*log10(abs(fftShifted)/max(abs(fftShifted)));
end

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

figure(1); plot(omega1khz, fftshiftNormdB(fft1));
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
figure(2); plot(omega1khz, fftshiftNormdB(fft1LPF)); title("lpf");

%% Quantization of low power elements
fft1Quant = fft1;
threshold = 0.00001 * max(abs(fft1));
lowIndex = abs(fft1) < threshold;
fft1Quant(lowIndex) = 0;

figure(3); plot(omega1khz, fftshiftNormdB(fft1Quant)); title("quantised");
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
prev = zeros(1,length(fft1));
errors = zeros(1,length(thresholds));
coeffDiff = zeros(1,length(thresholds));
numCoeffDeleted = zeros(1,length(thresholds));
for i = 1:length(thresholds)
    fft1Quant = fft1;
    threshold = thresholds(i) * max(abs(fft1));
    lowIndex = abs(fft1) < threshold;
    fft1Quant(lowIndex) = 0;

    % ifft & error
    seg1Out = ifft(fft1Quant);
    error = (seg1Out - seg1);
    % Ignore 1st/last 160 values as there are edge artefacts
    % av_error = mean(abs(error(160:end-160)));
    ssError = sum(error(160:end-160).^2);
    errors(i) = ssError;
    coeffDiff(i) = sum(lowIndex-prev);
    prev = lowIndex;
    numCoeffDeleted(i) = sum(lowIndex);
    % errors(floor(i/0.00005)+1) = av_error;
end

figure(5);
plot(thresholds*max(abs(fft1)), 20*log10(abs(errors)))
title("Normalised Error by Threshold Level")
xlabel("Lowest magnitude kept")
ylabel("Normalised Error (dB)")

figure(6);
plot(numCoeffDeleted, 20*log10(abs(errors)))
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
plot(omega1khz, fftshiftNormdB(fft1Quant));
title("threshold = 0.00062");
temp1 = fft1Quant;

%% Dynamic Thresholds
% Base thresholds around peaks
function [fftFull, dynThresholds] = dynMasking(fftSig, maskedMult)
    % Identify peaks in the original FFT (only in positive frequencies, 
    % just replicate it for negative frequencies later)
    fftMid = floor(length(fftSig)/2);
    fftHalf = fftSig(1:fftMid+1);
    fftdBHalf = 20*log10(abs(fftHalf)/max(abs(fftHalf))+1e-12);
    [peaks, locs] = findpeaks(fftdBHalf, "MinPeakProminence",10);
    
    % Define rate of decay of masking constants
    upwardsDecay = 0.5; % Masking decays slower upwards
    downwardsDecay = 2;
    
    % Calculate dynamic thresholds based on peak and their locations
    dynThresholds = -inf(size(fftdBHalf));
    for i = 1:length(locs)
        baseLoc = locs(i);
        
        for bin = 1:fftMid+1
            displacement = bin - baseLoc;
    
            % Calculate decay based on distance. Higher distance = less 
            % masking
            if (displacement > 0)
                decay = displacement * upwardsDecay;
            else 
                decay = -displacement * downwardsDecay;
            end
            
            % Use the least overall decay/highest threshold (i.e. apply the
            % most masking)
            dynThresholds(bin) = max(dynThresholds(bin), peaks(i) - decay);
        end
    end
    
    %% Apply dynamic thresholds to the FFT
    % Still simply cut off the frequency, but at least it is according to 
    % nearby frequencies
    % fft1dBHalf(fft1dBHalf < dynThresholds) = -Inf;
    fftHalf(fftdBHalf < dynThresholds) = fftHalf(fftdBHalf < dynThresholds).*abs(maskedMult);
    
    % Flip and add in the negative frequencies
    % fft1dBFull = [fft1dBHalf, flip(fft1dBHalf(2:end-1))];
    fftFull = [fftHalf, conj(flip(fftHalf(2:end-1)))];
end

%% Smaller windows
% Split signal into smaller overlapping blocks and apply the thresholds 
% for each block
winSize = 2048; % Size of block
overlap = 0.5; % Amount of overlap
maskedMod  = 0.5; % value to reduce masked frequencies (multiply)

function [sigOut] = windowSignal(signal, fftLength, winSize, overlap, maskedMod)
    winStep = winSize * overlap; % Distance between block centers
    hWindow = transpose(hann(winSize, "periodic")); % hann window
    winSums = zeros(size(signal)); % sum of hann windows to normalise
    sigOut = zeros(size(signal)); % output
    
    for step = 1:winStep:(fftLength - winSize + 1)
        block = signal(step:(step + winSize - 1)) .* hWindow;
    
        % Process the next block of the signal
        fftBlock = fft(block);
    
        % Apply dynamic masking to the FFT of the current block
        [fftBlock, ~] = dynMasking(fftBlock, maskedMod);
    
        % Add filtered data together
        sigOut(step:(step + winSize - 1)) = sigOut(step:(step + winSize - 1)) + ifft(fftBlock) .* hWindow;
    
        % Gather normaliser values
        winSums(step:(step + winSize - 1)) = winSums(step:(step + winSize - 1)) + hWindow .^2;
    end
    
    sigOut = sigOut ./ (winSums + 1e-12);
end

seg1Out = windowSignal(seg1, length(fft1), winSize, overlap, maskedMod);
sound(seg1Out(1:10*fs), fs);

%%
figure(9);
% plot(omega1khz, fftshiftNormdB(10.^(fft1dBFull./20)));
% title("dynamic threshold");
% 
% seg1Out = ifft(10.^(fft1dBFull./20));
% sound(seg1Out(1:7*fs), fs);
plot(omega1khz, fftshiftNormdB(fft(seg1Out)));
title("dynamic threshold");

% seg1Out = ifft(fft1Full);
% sound(seg1Out(1:7*fs), fs);
