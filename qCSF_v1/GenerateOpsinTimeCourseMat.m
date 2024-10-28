clear all; close all; clc;

addpath(genpath('../../UWToolbox'));
addpath(genpath('../utilities/'));
addpath(genpath('../model/'));
addpath(genpath("../"));
dataDir = ["../Data/LoadMats/"];
if ~exist(dataDir)
    mkdir(dataDir)
end

sz_deg = [8,8]; % size in degrees - this has to match the experiment parameters
display = struct('resolution', [1920 1080], ...
    'width', 70, 'dist', 206, ...
    'screen', 0, ...
    'fRate', 120, ...
    'pxPerDeg', 0, ...
    'anglePerPix', 0, ...
    'gray', 0);

matSz_pix = angle2pix(display,sz_deg(2));
sz = [matSz_pix, matSz_pix];   % size (pixels)
ar = sz(1)/sz(2);   % aspect ratio
ang = -45;           % orientation of Gabor
sigma = 1;         % size of Gaussian for the Gabor
contrast = 1;      
dt = 20;           % time step (msec) should be frame rate
dur = 800;            % duration (milliseconds)
pad = 1000;

Fixed_TF = [2, 4, 8, 12, 16]; %16*(1/2).^(4:-1:0);
Variable_SF = 10.^linspace(log10(.25),log10(16),12)';

Fixed_SF = exp(linspace(log(0.5), log(16),5));
Variable_TF = 10.^linspace(log10(3),log10(18),12)';

TFs = [Fixed_TF'; Variable_TF]; 
nC = 1024;  %number of contrasts in bank
cList = linspace(-1,1,nC);
t = 0:1/dt:(dur);
t_combined = 0:1/dt:(dur+pad);

lum = 0.5; opsin = ["ChRmine"]; % 'ChR2', 'ReaChR', 'ChrimsonR', 'CsChrimson', 'bReaChES', 'ChRmine'

% for i = 1:length(Variable_SF)
%     sf = Variable_SF(i);
%     % make the static image to be modulated (Gabor)
%     [x,y] = meshgrid(linspace(-1,1,sz(2)),linspace(-ar,ar,sz(1)));
%     ramp = cos(ang*pi/180)*x + sin(ang*pi/180)*y;
% %     Gauss = exp(-(x.^2+y.^2)/(2*sigma^2));
%     img  = contrast*cos(2*pi*sf*ramp);%.*Gauss;

    for j = 1:length(TFs)
        tf = TFs(j);
%     tf = Fixed_TF(j);        dataDir = ["../Data/LoadMats/"];
        if ~exist(dataDir)
            mkdir(dataDir)
        end
        if ang < 0
            nameang = 180-ang;
        else
            nameang = ang;
        end

        name = [dataDir + "tf_" + tf + ".mat"];
        % name = [dataDir + "temp2.mat"];
        if 1%~exist(name)
        % Generate the banks of contrast and opsin modulated time courses
        bank = [ones(dt*pad,1024)/2; ((1+sin(2*pi*tf*t/1000)'*cList)/2)];
        y = zeros(size(bank));
        for k = 1:length(cList)
            V = -60; lambda = 590; 
            LuxScaleFactor = (lum)*luxtoirradiance(lambda*1e-3);
            opsin_model = opsin_photocurrent();
            y(:,k) = opsin_model.get_opsin_current(opsin, V, lambda, LuxScaleFactor*squeeze(bank(:,k)));
        end

%         y = bank; % identity matrix for sanity check
        if tf == min(Fixed_TF)
            opto = p2p_opto(opsin);
            [offset, scaleFac] = opto.get_scalefactor(y(:,1024), tf, lum);
        elseif ~exist('offset', 'var')
            opto = p2p_opto(opsin);
            [offset, scaleFac] = opto.get_scalefactor(y(:,1024), 2, 0.5);
        end
        
        
        y = (y - offset)*scaleFac;
        y = 0.5*(y+1);

        figure
        plot(t,(bank(dt*pad+1:end,1024)));
        hold on;
        plot(t,(y(dt*pad+1:end,1024)));
        grid on;

        count=0;
        for g = dt*pad+1:166:length(t_combined)
            temp = y(g,:);
            count = count+1;
            SkipFramesBank(count,:) = temp;

        end
        
        % max(max(SkipFramesBank(:,1024)))
        save(name, 'SkipFramesBank');
        else
            load(name, 'SkipFramesBank');
        end
        
%         figure
%         n = linspace(0,800, size(SkipFramesBank,1));
%         plot(n,(SkipFramesBank(:,1024)))
% 
%         
%         tSeries = zeros(size(SkipFramesBank,1),prod(sz));
%         tic
%         img = img(:);
%         for g=1:(nC-1)
%             id = img >= cList(g) & img<=cList(g+1);
%             tSeries(:,id) = repmat(SkipFramesBank(:,g),1,sum(id));
%         end
%         toc
%         
% 
%         
% %         show the movie
% 
%         figure(3)
%         clf
%         curFrame = (reshape(tSeries(1,:),sz)+1)*128;
%         h = image(curFrame);
%         colormap(gray(256))
%         axis equal
%         axis off
%         truesize
% 
%         tic
%         for g=1:size(tSeries,1)
%             curFrame = (reshape(tSeries(g,:),sz)+1)*128;
%             set(h,'CData',curFrame)
%             drawnow
%             while(toc<t(g))
%                 % do nothing
%             end
%         end
%         toc
%             
    end    

% end