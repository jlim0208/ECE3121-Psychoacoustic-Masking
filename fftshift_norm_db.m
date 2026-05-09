function [output] = fftshift_norm_db(fftSig)
% FFTHIFT_NORM_DB Shift fft to 0-center, normalise and use dB scale
% Inputs:
%   fftSig      - Precomputed FFT of the full signal (vector)
% 
% Outputs:
%   output      - Normalised magnitude spectrum in dB with zero-frequency centered.
%                 Values are relative to the maximum magnitude (0 dB).
%
% Behavior:
%   The function takes a precomputed FFT vector (fftSig), computes its magnitude,
%   shifts the zero-frequency component to the center using fftshift, normalises
%   the magnitudes by the maximum magnitude to produce values in [0, 1], and
%   converts to decibels using 20*log10. If the input is all zeros, the output
%   will be -Inf due to log of zero.
%
% Notes:
%   - Function documentation generated using copilot and modified
%
% Example:
%   % Compute FFT of a signal and display centred, normalised dB spectrum
%   x = sin(2*pi*50*(0:1/1000:1-1/1000));
%   X = fft(x);
%   dbSpec = fftshift_norm_db(X);
%   f = (-floor(length(X)/2):ceil(length(X)/2)-1) * (1000/length(X));
%   plot(f, dbSpec); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
    fftShifted = fftshift(abs(fftSig));
    output = 20*log10(abs(fftShifted)/max(abs(fftShifted)));
end