function scale_factor = luxtoirradiance(lambda)
% lambda should be input in micrometers

% V(lambda) = 1.019*exp(-285.4* (lambda-0.559)^2)
% where, lambda is the wavelength of the monochromatic light - according to
% CIE and ISO standards
% Ref: https://www.thorlabs.de/catalogPages/506.pdf
% There is 683 lumens per watt for 555 nm photopic light

V = 1.019*exp(-285.4* (lambda-0.559)^2);
scale_factor = 1/(V*683);
end