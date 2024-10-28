classdef p2p_opto

    properties
        opsin = "ChRmine";
    end
methods (Access = public)

    function this = p2p_opto(opsin)
        this.opsin = opsin; % ["ChR2", "ReaChR", "ChrimsonR", "CsChrimson", "bReaChES", "ChRmine"];
    end

    function [offset, scaleFac] = get_scalefactor(this, y, min_TF, luminance)
%         V = -60; lambda = 590;
        pad = 1000; % Specify the amount of time to pad with ones before and after stimulus in ms 
        dt = 20;
%         time = 0:1/dt:(1000 + pad);
% 
%         LuxScaleFactor = (luminance)*luxtoirradiance(lambda*1e-3);
%         Irr = LuxScaleFactor*[ones(1, dt*pad)/2, (1 + sin(2*pi*min_TF*(time((dt*pad+1):end))/1000))/2]';
%         opsin_model = opsin_photocurrent();
%         y = opsin_model.get_opsin_current(this.opsin, V, lambda, Irr);


        y0 = y(dt*pad+1);
        a = max(y(dt*pad+1:end));
        b = min(y(dt*pad+1:end));
        scaleFac = 1/(max(a-y0, y0-b)); % This scales y such that it varies between 0 and 1
        offset = y0; % this offsets y such that it begins at y0
    end
end
end