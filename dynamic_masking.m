function [fftOutFull, lowIndex] = dynamic_masking(fftSig, maskedMult, decays)
% DYNAMIC_MASKING Apply masking around peaks
% Inputs:
%   fftSig     - Precomputed FFT of the full signal (vector)
%   maskedMult - Scalar multiplier applied to masked bins (typically <1)
%   decays     - Two-element vector specifying asymmetric decay rates (in 
%                dB/bin) for masking thresholds away from each detected 
%                peak. Follow format of: decay for frequencies below, then 
%                above. Larger values produce faster reduction of the 
%                masking threshold with distance. Defaults to [2, 0.5].
%
% Outputs:
%   fftOutFull - FFT after applying dynamic masking (same size as fftSig)
%   lowIndex   - Logical vector indicating bins that were masked (true)
%
% Behavior:
%   The function identifies prominent peaks in the positive-frequency
%   half of the FFT, builds per-bin dynamic masking thresholds that
%   decay away from each peak (asymmetric decay rates for upwards and
%   downwards frequency directions), and attenuates FFT bins whose
%   magnitudes (in dB relative to the largest bin) fall below those
%   thresholds. Masking is applied symmetrically to the negative
%   frequency bins by mirroring the positive-half results.
%
% Notes:
%   - The input fftSig is expected to be a full FFT vector (length N).
%   - maskedMult should be a scalar. If not provided, it defaults to 0.5.
%
% Validate inputs
if nargin < 1
    error('dynamic_masking requires at least one input: fftSig');
end
if nargin < 2 || isempty(maskedMult)
    maskedMult = 0.5; % default to half attenuation of masked bins
end
if ~isvector(fftSig)
    error('fftSig must be a vector.');
end
if ~isscalar(maskedMult) || ~isnumeric(maskedMult)
    error('maskedMult must be a numeric scalar.');
end
if nargin < 3 || isempty(decays)
    decays = [2, 0.5];
end
if iscell(decays)
    decays = cell2mat(decays);
end

% Base thresholds around peaks which effect decays as you get further away
% from them
    % Identify peaks in the original FFT (only in positive frequencies, 
    % just replicate it for negative frequencies later)
    fftMid = floor(length(fftSig)/2);
    fftHalf = fftSig(1:fftMid+1);
    fftdBHalf = 20*log10(abs(fftHalf)/max(abs(fftHalf))+1e-12);
    [peaks, locs] = findpeaks(fftdBHalf, "MinPeakProminence",10);
    
    % Define rate of decay of masking constants
    upwardsDecay = decays(2); % Masking decays slower upwards
    downwardsDecay = decays(1);
    
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
    lowIndexHalf = fftdBHalf < dynThresholds;
    fftHalf(lowIndexHalf) = fftHalf(lowIndexHalf).*abs(maskedMult);
    
    % Flip and add in the negative frequencies
    fftOutFull = [fftHalf, conj(flip(fftHalf(2:end-1)))];
    lowIndex = [lowIndexHalf, flip(lowIndexHalf(2:end-1))];
end