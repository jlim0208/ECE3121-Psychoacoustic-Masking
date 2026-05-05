%% Test multiple thresholds
% Test different threshold levels using static_masking
function [errors, coeffDiff, numCoeffDel] = test_static_thresholds(fft, thresholds, ogSignal)
    % Preallocate
    prev = zeros(1,length(fft)); % Save previous coefficients deleted to 
                                  % calculate difference
    errors      = zeros(1, length(thresholds));
    coeffDiff   = zeros(1, length(thresholds));
    numCoeffDel = zeros(1, length(thresholds));
    
    % Main testing loop
    for i = 1:length(thresholds)
        [fftQuant, lowIndex] = static_masking(fft, thresholds(i) * max(abs(fft)));
    
        % ifft & error calculations
        sigOut = ifft(fftQuant);
        error = (sigOut - ogSignal);
        ssError = sum(error(160:end-160).^2); % Ignore 1st/last 160 values 
                                              % as there are edge artefacts
        errors(i) = ssError;
    
        coeffDiff(i) = sum(lowIndex-prev);
        prev = lowIndex; 
        numCoeffDel(i) = sum(lowIndex);
    end
end