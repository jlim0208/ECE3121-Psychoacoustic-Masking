%% Smaller windows
% Split signal into smaller overlapping blocks and apply dynamic thresholds
% for each block
function [sigOut] = window_signal(signal, fftLength, winSize, overlap, maskedMod)
    winStep = winSize * overlap; % Distance between block centers
    hWindow = transpose(hann(winSize, "periodic")); % hann window
    winSums = zeros(size(signal)); % sum of hann windows to normalise
    sigOut = zeros(size(signal)); % output
    
    for step = 1:winStep:(fftLength - winSize + 1)
        block = signal(step:(step + winSize - 1)) .* hWindow;
    
        % Process the next block of the signal
        fftBlock = fft(block);
    
        % Apply dynamic masking to the FFT of the current block
        [fftBlock, ~] = dynamic_masking(fftBlock, maskedMod);
    
        % Add filtered data together
        sigOut(step:(step + winSize - 1)) = sigOut(step:(step + winSize - 1)) + ifft(fftBlock) .* hWindow;
    
        % Gather normaliser values
        winSums(step:(step + winSize - 1)) = winSums(step:(step + winSize - 1)) + hWindow .^2;
    end
    
    sigOut = sigOut ./ (winSums + 1e-12);
end