%% Kelly_1979_stuff.m
%
% Generates fake psychophysical 2AFC data from the parameterized Kelly TCSF
% and reconstructs the surface from the fake data.
%
% gmb 4/5/2024
%
% See also 'KellySurface.m', 'fitKelly', 'weibull.m'
clear all; clc

cond  = {'baseline' 'opto' 'eye'};
for i = 1:length(cond)

    itercond = cell2mat(cond(i));
exptData = load(strcat("../analysis_data/NS_L_GED_1968_", itercond, ".mat"));
end
%%
% grid of spatial and temporal frequencies for plotting the smooth surfaces
[plotsf,plottf] = meshgrid(2.^linspace(0.5,5,101));

% initial (default) D.H. Kelly parameters
p.param1 = 6.1;
p.param2 = 7.3;
p.param3 = 45.9;
p.b = 2;  % slope of psychometric function

% find and plot the sensitivity surface
Sinit = calculateKelly(p,plotsf,plottf);

figure(1)
clf
surf(log(plotsf(1,:)),log(plottf(:,1)),Sinit);
set(gca,'Zlim',[0,2.5])
logx2raw
logy2raw
xlabel('Spatial Frequency')
ylabel('Temporal Frequency')


%%
% Generate fake psychophysical data:
% stimulus parameters sf, tf, contrast
% data parameter resp (0 or 1) 
%
% replace 'sf, tf, contrast and resp' with your data!

%  Stimulus parameters.  
% fList = 2.^[0.5:2:25];  % list of spatial and temporal frequencies
% cList = 2.^linspace(-7,-2,8); % list of contrasts
% nReps = 10;  % trials per sf,tf, and contrast

% dims = [length(fList),length(fList),length(cList),nReps];
% disp(sprintf('%d trials.',prod(dims)))

% create column vectors of stimulus parameters
% [sf,tf,contrast] = meshgrid(fList,fList,cList);

% sf = repmat(sf,[1,1,1,nReps]); sf = sf(:);
% tf = repmat(tf,[1,1,1,nReps]); tf = tf(:);
% contrast = repmat(contrast,[1,1,1,nReps]); contrast = contrast(:);

iterData = strcat(exptData.itercond);
sf =  exptData.opto(:,2);
tf =  exptData.opto(:,1);
contrast =  exptData.opto(:,3)./100;
resp = exptData.opto(:,4);

% Find the sensitivity for each stimulus
S = calculateKelly(p,sf,tf);

% convert from log10 sensitivity to threshold contrast.  Save them as
% threshold values for the Weibull function
p.t = .1.^S;    % or 1./S.^10

% Feed these thresholds into the Weibull function to get the expected
% P(correct) for each trial
% prob = weibull(p,contrast(:));

% fake psychophysical binomial responses based on biased prob:
% resp = floor(rand(size(prob))+prob);


%% 
% For fun, show the expected probability correct as a surface for a given
% contrast value

% tmpProb=  reshape(prob,dims);
% cNum = 3;  % choose a contrast in the list to display
% figure(2)
% clf
% surf(log(fList),log(fList),squeeze(tmpProb(:,:,cNum,1)));
% logx2raw
% logy2raw
% xlabel('Spatial Frequency')
% ylabel('Temporal Frequency')
% zlabel('P(correct)')
% title(sprintf('Contrast %5.2f%%',100*cList(cNum)))

%%
% Fit the Kelly surface to the psychophysical data using -log likelihood.

% Here's the initial error
errInit = fitKelly(p,sf(:),tf(:),contrast(:),resp);

figure(1)
title(sprintf('initial: p1 = %5.2f, p2= %5.2f, p3 = %5.2f, err= %5.2f',p.param1,p.param2,p.param3,errInit))

% Use fmins to find the best-fitting parameters
freeParams = {'param1','param2','param3'};
pBest = fit('fitKelly',p,freeParams,sf(:),tf(:),contrast(:),resp);
errBest = fitKelly(pBest,sf(:),tf(:),contrast(:),resp)

%%
% Plot the best-fitting surface.  Hopefully it looks like figure 1

Sbest = calculateKelly(pBest,plotsf,plottf);

figure(3)
clf
surf(log(plotsf(1,:)),log(plottf(:,1)),Sbest);
logx2raw
logy2raw
xlabel('Spatial Frequency')
ylabel('Temporal Frequency')
set(gca,'Zlim',[0,2.5])
title(sprintf('fit: p1 = %5.2f, p2= %5.2f, p3 = %5.2f, err= %5.2f',pBest.param1,pBest.param2,pBest.param3,errBest))

%%
% More plotting fun: plot the initial and best-fitting surfaces in the same figure.
figure(4)
clf
hold on

surf(log(plotsf(1,:)),log(plottf(:,1)),Sinit, 'FaceColor','g', 'FaceAlpha',0.2,'EdgeColor','none');
surf(log(plotsf(1,:)),log(plottf(:,1)),Sbest, 'FaceColor','b', 'FaceAlpha',.2,'EdgeColor','none');

set(gca,'Zlim',[0,2.5])
logx2raw
logy2raw
xlabel('Spatial Frequency')
ylabel('Temporal Frequency')

hold off
view(25,25)
legend({'init','fit'})