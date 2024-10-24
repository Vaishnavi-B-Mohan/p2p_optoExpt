close all; clear all; clc;

%% FitPsycho2Data.m
%
% explores fitting Kelly's surface to quickCSF data

% Col-1: TFs; Col-2: SFs; Col-3: Sensitivity values; Col-4: Responses.
% Intial Weibull parameters
initp.t = .15;
initp.b = 1;

ticks = [1,2,4,8,16,32];

% Grid in tf sf space
ntf = 4;
nsf = 4;
tfList = linspace(.4,3.8,ntf+1);
sfList = linspace(-1.5,3.25,nsf+1);
[tfGrid,sfGrid] = meshgrid(tfList,sfList);

% for plotting
plottf = exp((tfList(1:end-1)+tfList(2:end))/2);
plotsf = exp((sfList(1:end-1)+sfList(2:end))/2);

% get list of files
subjects = dir('../analysis_data/NS*.mat');

colList = {[1,0,0],[0,1,0],[0,0,1]};

S = zeros(ntf,nsf,length(subjects),3);

for subNum = 1:length(subjects)
    load(['../analysis_data/',subjects(subNum).name])
%     figure(subNum)
%     clf
%     hold on
    for expNum = 1:2

        % bin data by grid of tf and sf and fit psychometric functions for
        % each cell in the grid.

        for i=1:ntf
            for j=1:nsf
                id = log(exptData(expNum).data(:,1)) > tfList(i) & ...
                    log(exptData(expNum).data(:,1)) <= tfList(i+1) & ...
                    log(exptData(expNum).data(:,2)) > sfList(j) & ...
                    log(exptData(expNum).data(:,2)) <= sfList(j+1);

                clear results
                results.intensity = .1.^exptData(expNum).data(id,3);  %YES!
                results.response = logical(exptData(expNum).data(id,4));

                pBest = fitcon('fitWeibull',initp,{'t>0'},results);
                S(i,j,subNum,expNum) = -log10(pBest.t);
            end
        end
%        surf(log(plotsf),log(plottf),S(:,:,subNum,expNum), 'FaceColor',colList{expNum}, 'FaceAlpha',0.2,'EdgeColor','none');

    end
%     logx2raw
%     logy2raw
%     xlabel('Spatial Frequency')
%     ylabel('Temporal Frequency')
%     set(gca,'Zlim',[0,2.5])
% 
%     view(25,25)
xlabel('Spatial Frequency')
size(S)
end

%%
% show the sensitivity surfaces averaged across subjects
Smean = squeeze(mean(S,3));

figure(length(subjects)+1)
clf
hold on
for expNum = 1:2
    surf(log(plotsf),log(plottf),Smean(:,:,expNum), 'FaceColor',colList{expNum}, 'FaceAlpha',0.2,'EdgeColor','none');
end
set(gca,'XTick',log(ticks))
set(gca,'YTick',log(ticks))
set(gca,'XTickLabel',num2str(ticks'))
set(gca,'YTickLabel',num2str(ticks'))

legend()
% logx2raw
% logy2raw
xlabel('Spatial Frequency')
ylabel('Temporal Frequency')
legend({'control','opto','eye'})
view(25,25)

% rotate!
axis vis3d
grid
for i=1:10
    for ang = 1:.25:360
        view(ang,25)
       
        drawnow
    end
end

