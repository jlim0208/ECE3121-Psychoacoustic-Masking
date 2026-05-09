function [fftOut, lowIndex] = window_signal(signal, fftSig, winSize, overlap, maskParam, masking_func)
% WINDOW_SIGNAL Apply masking by windowing
%
% Inputs:
%   signal       - Time-domain signal vector to be processed (column or row).
%   fftSig       - Precomputed FFT of the full signal (vector). Its length
%                  determines processing range.
%   winSize      - Size (in samples) of each window/block (positive integer).
%   overlap      - Fractional overlap between windows expressed as a scalar
%                  step multiplier (0 < overlap <= 1). winStep is computed
%                  as winSize * overlap (should produce integer step).
%   maskParam    - Additional parameter passed to masking_func (e.g. mask
%                  threshold or attenuation factor).
%   masking_func - Function handle: [fftFiltered, lowIndex] = masking_func(fft, maskedMod)
%                  that applies masking to a block FFT and returns the
%                  masked FFT and a logical or numeric vector indicating
%                  which frequency bins were quantised.
%
% Outputs:
%   fftOut   - Reconstructed frequency-domain signal after block-wise 
%              masking and overlap-add (same size as input signal).
%   lowIndex - Row vector containing, for each processed block, the count
%              of quantised/low-resolution bins returned by masking_func.
%
% Behavior:
%   Split signal into smaller overlapping blocks and apply masking for each 
%   block separately, then combine them using hann windows.
%
% Notes:
%   - The function uses a periodic Hann window for analysis/synthesis and
%     normalises the overlap-add using the sum of squared windows.
%   - It is assumed that winSize and winStep (winSize * overlap) are
%     chosen so that window positions align within the length of fftSig.
%   - The small constant added during division avoids numerical issues.
%   - Function documentation generated using copilot and modified
%
% Example:
%   X = fft(x);
%   winSize = 2048;
%   overlap = 0.5;
%   maskedMod  = 0.5;
%   [sigOut, lowIndex] = window_signal(x, X, winSize, overlap, maskedMod, @dynamic_masking)

    fftLength = length(fftSig);
    winStep = winSize * overlap; % Distance between block centers
    hWindow = transpose(hann(winSize, "periodic")); % hann window
    winSums = zeros(size(signal)); % sum of hann windows to normalise
    sigOut = zeros(size(signal)); % output
    numSteps = floor((fftLength - winSize)/winStep) + 1;
    lowIndex = zeros(1, numSteps); % Each element is the number of 
                                   % quantised elements in each block
    count = 1;    
    
    for step = 1:winStep:(fftLength - winSize + 1)
        block = signal(step:(step + winSize - 1)) .* hWindow;
    
        % Process the next block of the signal
        fftBlock = fft(block);
    
        % Apply masking to the FFT of the current block
        [fftBlock, lowIndexBlock] = masking_func(fftBlock, maskParam);
    
        % Add filtered data together
        sigOut(step:(step + winSize - 1)) = sigOut(step:(step + winSize - 1)) + ifft(fftBlock) .* hWindow;
    
        % Gather normaliser values
        winSums(step:(step + winSize - 1)) = winSums(step:(step + winSize - 1)) + hWindow .^2;

        % Sum number of quantised bins
        lowIndex(count) = sum(lowIndexBlock);
        count = count + 1;
    end
    
    sigOut = sigOut ./ (winSums + 1e-12);
    fftOut = fft(sigOut);
end