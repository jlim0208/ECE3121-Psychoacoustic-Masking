function [errors, coeffDiff, numCoeffDel] = test_thresholds(fftSig, thresholds, ogSignal, masking_func)
% TEST_THRESHOLDS test multiple variables and accumulate error in masking
% Inputs:
%   fftSig      - Precomputed FFT of the full signal (vector)
%   thresholds  - Vector of threshold values to test with masking_func
%   ogSignal    - Original time-domain signal corresponding to fft (vector)
%   masking_func- Function handle: [fftQuant, lowIndex] = masking_func(fft, threshold)
%                 where fftQuant is the thresholded FFT and lowIndex is a logical
%                 or numeric mask (1 for removed/low coefficients, 0 otherwise)
%
% Outputs:
%   errors      - Vector of sum-squared reconstruction errors for each threshold
%   coeffDiff   - Vector of differences in number of removed coefficients
%                 compared to the previous threshold (first entry is removed)
%   numCoeffDel - Vector of total number of removed coefficients for each threshold
%
% Notes:
%   - The function compares the ifft of the masked FFT with the original signal and
%     computes the sum of squared errors, ignoring the first and last 160 samples
%     to avoid edge artifacts.
%   - masking_func is expected to return consistent-sized outputs and a mask that
%     can be summed to count removed coefficients.
%   - Function documentation generated using copilot and modified

    % Preallocate
    prev = zeros(size(fftSig)); % Save previous coefficients deleted to 
                             % calculate difference
    errors      = zeros(size(thresholds));
    coeffDiff   = zeros(size(thresholds));
    numCoeffDel = zeros(size(thresholds));
    
    % Main testing loop
    for i = 1:length(thresholds)
        [fftQuant, lowIndex] = masking_func(fftSig, thresholds(i));
    
        % ifft & error calculations
        sigOut = ifft(fftQuant);
        error = (sigOut - ogSignal);
        ssError = sum(error(160:end-160).^2); % Ignore 1st/last 160 values 
                                              % as there are edge artefacts
        errors(i) = ssError;
    
        coeffDiff(i) = sum(lowIndex)-sum(prev);
        prev = lowIndex; 
        numCoeffDel(i) = sum(lowIndex);
    end
end