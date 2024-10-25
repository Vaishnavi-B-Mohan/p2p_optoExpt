function [phaseshifted_grating] = changephase_fft(grating, f_ang_unshifted)
f = fft2(grating);
% scaleFac = amplitude_unshifted/max(max(abs(f)));
% new_fft = pol2cart(f_ang_unshifted, scaleFac.*abs(f));
new_fft = pol2cart(f_ang_unshifted, abs(f));
phaseshifted_grating = ifft2(new_fft);
end