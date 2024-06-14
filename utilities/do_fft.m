function [f_ang_unshifted] = do_fft(grating)
f=fft2(grating);
% Shift FT from corners to central part.
% f_mag =log(abs(fftshift(f)));
% f_ang = (angle(fftshift(f)));
f_ang_unshifted = (angle(f));
% figure
% subplot(2,1,1)
% imshow(f_mag);
% title('FFT magnitude spectrum')
% subplot(2,1,2)
% imshow(f_ang)
% title('FFT phase spectrum')
end