% This is a function to convert pixels to visual angle given the dimensions
% of the display and it's resolution
% Inputs: 
% 1) dim: Dimensions of the display and the viewing distance (all in the same metric, could be inches. cm etc)
%         [width, height, depth]
% 2) rn: Flag to indicate whether the visual angle must be returned in
% radians (=1) or degrees (=0). Default is in radians (=1)

% Output(s):
% pix_angle: Vector specifying the visual span of the display in visual
% angle [horizontal, vertical]

function [pix_angle] = PixelToVisAngle(dim, rn)
    if (~exist('rn', 'var'))
        rn = 1;
    end
    pix_angle = zeros(2,1);
    if (rn == 1)
        pix_angle(1) = 2*atan(dim(1)/(2*dim(3)));
        pix_angle(2) = 2*atan(dim(2)/(2*dim(3)));
    else if (rn == 0)
        pix_angle(1) = (180/pi)*2*atan(dim(1)/(2*dim(3))); 
        pix_angle(2) = (180/pi)*2*atan(dim(2)/(2*dim(3)));
    else
        print("Parameter rn only takes values 0 or 1\n");
        pix_angle = [];
    end
end