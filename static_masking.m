function [fftOut, lowIndex] = static_masking(fftSig, threshold)
% STATIC_MASKING Removes all signals lower than a given threshold
% Inputs:
%   fftSig    - Precomputed FFT of the full signal (vector)
%   threshold - Scalar magnitude threshold. Frequency components with magnitude
%               less than this value are considered noise and will be zeroed.
% 
% Outputs:
%   [fftOut, lowIndex] = static_masking(fftSig, threshold) returns the
%   masked FFT (fftOut) and a logical mask (lowIndex) indicating which
%   frequency bins were below the threshold.
%
%   Example:
%     X = fft(x);
%     [Xmask, idx] = static_masking(X, 1e-3);
%
%   The function operates element-wise on fftSig and supports complex inputs.
    fftOut = fftSig;
    lowIndex = abs(fftSig) < threshold;
    fftOut(lowIndex) = 0;
end