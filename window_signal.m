%% Smaller windows
% Split signal into smaller overlapping blocks and apply masking for each 
% block separately, then combine them using hann windows
function [sigOut, lowIndex] = window_signal(signal, fftSig, winSize, overlap, maskedMod, masking_func)
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
        [fftBlock, lowIndexBlock] = masking_func(fftBlock, maskedMod);
    
        % Add filtered data together
        sigOut(step:(step + winSize - 1)) = sigOut(step:(step + winSize - 1)) + ifft(fftBlock) .* hWindow;
    
        % Gather normaliser values
        winSums(step:(step + winSize - 1)) = winSums(step:(step + winSize - 1)) + hWindow .^2;

        % Sum number of quantised bins
        lowIndex(count) = sum(lowIndexBlock);
        count = count + 1;
    end
    
    sigOut = sigOut ./ (winSums + 1e-12);
end