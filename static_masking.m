%% Apply a flat mask
% Removes all signals lower than a given threshold
function [fftOut, lowIndex] = static_masking(fftSig, threshold)
    fftOut = fftSig;
    lowIndex = abs(fftSig) < threshold;
    fftOut(lowIndex) = 0;
end