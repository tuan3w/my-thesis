function x = snr(s,r)
% SNR calculates the signal-to-noise ratio between the signal s and its
% representation r defined as:
% SNR = 20*log10(||s||_2/||s-r||_2)
err = s-r;
x = mag2db(norm(s(:))/norm(err(:)));