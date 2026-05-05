%% Shift fft to 0-center, normalise and use dB scale
function [output] = fftshift_norm_db(fft)
    fftShifted = fftshift(abs(fft));
    output = 20*log10(abs(fftShifted)/max(abs(fftShifted)));
end