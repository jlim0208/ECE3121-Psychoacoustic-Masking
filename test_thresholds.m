%% Test multiple thresholds
% Test different threshold levels using static_masking
function [errors, coeffDiff, numCoeffDel] = test_thresholds(fft, thresholds, ogSignal, masking_func)
    % Preallocate
    prev = zeros(size(fft)); % Save previous coefficients deleted to 
                             % calculate difference
    errors      = zeros(size(thresholds));
    coeffDiff   = zeros(size(thresholds));
    numCoeffDel = zeros(size(thresholds));
    
    % Main testing loop
    for i = 1:length(thresholds)
        [fftQuant, lowIndex] = masking_func(fft, thresholds(i));
    
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