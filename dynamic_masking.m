%% Dynamic Thresholds
% Base thresholds around peaks which effect decays as you get further away
% from them
function [fftOutFull, dynThresholds] = dynamic_masking(fftSig, maskedMult)
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
    fftOutFull = [fftHalf, conj(flip(fftHalf(2:end-1)))];
end