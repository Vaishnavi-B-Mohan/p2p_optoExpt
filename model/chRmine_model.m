% Theoretical models for response kinetics of various opsins based on the
% Fohlmeister and Miller model. Ref: https://iopscience.iop.org/article/10.1088/1741-2552/ac1175#references
% Written by Vaishnavi Mohan: November 2023
% Modified: August 2024
% Vision and Cognition Group, University of Washington

clear all; close all; clc;

%% Plots to visualize open and closed states
% plot_x = linspace(0,1,NumPts);
% figure(1)
% subplot(2,2,1)
% plot(plot_x, C1);
% ylabel('C1');
% subplot(2,2,2)
% plot(plot_x, O1);
% ylabel('O1');
% subplot(2,2,3)
% plot(plot_x, O2);
% ylabel('O2');
% subplot(2,2,4)
% plot(plot_x, C2);
% ylabel('C2');


%% Gating function parameters of ion-channels for RGCs
V = -60; % held at -60 mV

% Get ionic currents and visualize
FM_model = ionic_channel_fm();
I_ionic = get_ionic_current(FM_model, V);
dt = 0.001;
figure(2)
time = 0:dt:1;
plot(time, I_ionic);
grid on;
title('Ionic current: sum of all ion channels')

%% Opsin current
close all; clear all; clc;
V = -60;
% Get opsin current and visualize
opsins = ["ChR2", "ReaChR", "ChrimsonR", "CsChrimson", "bReaChES", "ChRmine"];
% opsins = ["ChRmine"];
lambda = 590; % wavelength of incident light in nm
pad = 1000; % Specify the amount of time to pad with ones before and after stimulus in ms 
dt = 50;
time = 0:1/dt:(1000 + 2*pad);
nt = length(time);
flist = [0.5, 2, 8, 16, 22, 30]; % temporal frequency of input in Hz
% flist = [1.875];
LuxScaleFactor = luxtoirradiance(lambda*1e-3);
Amp = 0.5*LuxScaleFactor; % in W/mm2
Iopsin_range = zeros(length(opsins), 1);
for i = 1:length(opsins)
    Iopsin_max = -inf; Iopsin_min = inf;
    for f = 1:length(flist)
        opsin_model = opsin_photocurrent();
        Irr = Amp*[ones(1, dt*pad)/2, (1 + cos(2*pi*flist(f)*(time((dt*pad+1):end-dt*pad))/1000))/2, ones(1, dt*pad)/2]';
        I_opsin = opsin_model.get_opsin_current(opsins(i), V, lambda, Irr);
        iter_max = max(I_opsin(dt*pad+1:end-dt*pad));
        iter_min = min(I_opsin(dt*pad+1:end-dt*pad));
        if iter_max > Iopsin_max
            Iopsin_max = iter_max;
        end
        if iter_min < Iopsin_min
            Iopsin_min = iter_min;
        end
        
        figure(i)
        subplot(3,2,f)
        plot(time(dt*pad-200:end-dt*pad+200), I_opsin(dt*pad-200:end-dt*pad+200));
%         plot(time, I_opsin);
%         x`ylim([Iopsin_min, Iopsin_max]);
        titlestr = ['Opsin: ' + opsins(i) + "  TF = " + flist(f)];
        title(titlestr);
        xlabel('Time in ms');
        grid on;
        xlabel('Time in ms');
        hold on;
        if i == length(opsins)
            figure(8)
            subplot(3,2,f)
            plot(time(dt*pad-200:end-dt*pad+200), Irr(dt*pad-200:end-dt*pad+200));
            title('Input irradiance  in W/mm2');
            grid on;
        end
%         Iopsin_max(i,f) = max(I_opsin(dt*pad+1:end-dt*pad));
%         Iopsin_min(i,f) = min(I_opsin(dt*pad+1:end-dt*pad));
    end
    Iopsin_range(i) = Iopsin_max  - Iopsin_min;
end

scale_factor = Iopsin_range./max(Iopsin_range);



% legend('ChR2', 'ReaChR', 'ChrimsonR', 'CsChrimson', 'bReaChES', 'ChRmine');
% Irr = [zeros(1, (nt-1)/10), Amp*square(2*pi*f*((time(1:(nt+1)/2))/1000)), zeros(1, (nt-1)/2 - (nt-1)/10)];
% opsins = ['ChR2', 'ReaChR', 'ChrimsonR', 'CsChrimson', 'bReaChES', 'ChRmine'];
% for f = 1:length(flist)
% Irr = Amp*[ones(1, dt*pad)/2, (1 + cos(2*pi*f*(time((dt*pad+1):end-dt*pad))/1000))/2, ones(1, dt*pad)/2]';
% figure(1)
% plot(time, Irr)
% grid on;
% Irr = square(f*time/1000);
