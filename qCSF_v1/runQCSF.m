function varargout = runQCSF(qcsf,varargin)
% 
% RUNQCSF   This program implements the critical trial-to-trial components of quick CSF method
% 
% For a qCSF application, the runQCSF program is used to call the subroutines that will 
% perform pre and post trial calculations, and plot the results of qCSF simulations and
% experiments. The calls to runQCSF are:
%       runQCSF(qcsf,'pretrial'); % Pre-trial analysis
%  
%       runQCSF(qcsf,'posttrial',response); % Update qCSF estimates
%                                                 based on the observer's response: 
%                                                 correct(response=1) 
%                                                 or incorrect(response=0) 
% 
%       runQCSF(qcsf,'plot experiment');     % Plots expt results 
%     
%
%
% For simulations: 
%
%      runQCSF(qcsf,'simulate'); % A weighted coin-flip is used to
%                                       simulate the observer's response
%  
%      runQCSF(qcsf,'plot simulation'); % Plots simulation results 
%    
% Version History
% Beta, November 5, 2009             fixed a bug that crashed experimental applications
%
% Beta, December 1, 2009             added two optimization parameters for users to change, if they wish:
%                              (1) - sample number for Monte Carlo sampling of the prior 
%                              and (2) the top nth-percentile (of informative stimuli) over which
%                              stimuli are chosen uniformly. The qCSF simulations in the JOV paper
%                              used qcsf.parameters.priorSamples = 50, and qcsf.parameters.optPercentile=10;(10th percentile)
%                              These parameters can be changed in setupQCSF.m 
% Version 1.0, June 30, 2010             


stringID = varargin{1};
 
if strcmp('PRETRIAL',upper(stringID)),
    
    [qcsf,nextFrequency,nextContrast]= preTrialAnalysis(qcsf);

    varargout{2}=nextFrequency;
    varargout{3}=nextContrast;
    
elseif strcmp('SIMULATE',upper(stringID)),
    
    gratingFrequency = varargin{2};
    gratingContrast = varargin{3};
    
    [qcsf,expectedCorrect,response] = simulateResponse(qcsf,gratingFrequency,gratingContrast);
    
    varargout{2}=response;
    varargout{3}=expectedCorrect;
   
elseif strcmp('POSTTRIAL',upper(stringID)),
    
    gratingFrequency = varargin{2};
    gratingContrast = varargin{3};
    response = varargin{4};
    
    qcsf=postTrialAnalysis(qcsf,gratingFrequency,gratingContrast,response);

elseif strcmp('PLOT SIMULATION',upper(stringID)), 
        
       qcsf=plotSimulation(qcsf);
       figure(gcf);
       
elseif strcmp('PLOT EXPERIMENT',upper(stringID)),  
          
        qcsf=plotExperiment(qcsf);
        figure(gcf);
       
end
       
varargout{1}=qcsf;
  
%%
function [qcsf,nextFrequency,nextContrast] = preTrialAnalysis(qcsf)

    priorSamples = qcsf.parameters.priorSamples;
    pTile = qcsf.parameters.optPercentile/100;
    
    if pTile < .01, pTile = .01;
    elseif pTile > 1, pTile = .99;
    end

    S = discretesamplev6(qcsf.data.prior,priorSamples);

    [sGain,sCenter,sWidth,sTrunc] = ind2sub(size(qcsf.data.prior),S);

    sCSF = [qcsf.parameters.gain(sGain)' qcsf.parameters.center(sCenter)' qcsf.parameters.width(sWidth)' qcsf.parameters.trunc(sTrunc)'] ;

    for nStimulus = 1:size(qcsf.stimuli.conditions,1),
        [psuccess,H]=findCSF_Probability(log10(qcsf.stimuli.conditions(nStimulus,1)),log10(qcsf.stimuli.conditions(nStimulus,2)),'column',...
            sCSF,qcsf.parameters.guessingRate,qcsf.parameters.lapseRate);
        I(nStimulus) = computeEntropy(mean(psuccess))-mean(H);
    end

    [maxI,indexI] = sort(I);
    indexI = indexI(end:-1:1);

    tmpI = randperm(ceil((pTile*length(maxI(:)))));
    indexI = indexI(tmpI);
    bestStimulus=indexI(1);

    nextFrequency = qcsf.stimuli.conditions(bestStimulus,1); 
    nextContrast = qcsf.stimuli.conditions(bestStimulus,2);

    qcsf.stimuli.nextFrequency=nextFrequency;
    qcsf.stimuli.nextContrast=nextContrast;
    qcsf.stimuli.nextEntropy = I;

function [qcsf, r, Response] = simulateResponse(qcsf,gratingFrequency,gratingContrast)
    
   r=findCSF_Probability(log10(gratingFrequency),log10(gratingContrast),'point',...
                    qcsf.simulation.trueCSF,qcsf.parameters.guessingRate,qcsf.parameters.lapseRate./qcsf.parameters.nAlternatives);

   Response = rand<r;

%%
function qcsf = postTrialAnalysis(qcsf,gratingFrequency,gratingContrast,response)

    prior=qcsf.data.prior;
    psuccess= findCSF_Probability(log10(gratingFrequency),log10(gratingContrast),'grid',...
                                qcsf,qcsf.parameters.guessingRate,qcsf.parameters.lapseRate./qcsf.parameters.nAlternatives);

    if response,
        pts_x = prior(1:end) *  psuccess(:);  
        pt2= (prior(1:end)').*psuccess(:)./pts_x;
    else
        pts_x = prior(1:end) *  (1-psuccess(:)); 
        pt2 = (prior(1:end)').*(1-psuccess(:))./(pts_x);
    end
   
    posterior = reshape(pt2,size(prior));
    qcsf.data.prior = posterior;
    
    mCSF= analyzePosterior(qcsf);
    qcsf.data.estCSF(qcsf.data.trial,:) = mCSF;
    qcsf.data.estSensitivity(qcsf.data.trial,:)=findQCSF(log10(qcsf.stimuli.frequency),mCSF(1),mCSF(2),mCSF(3),mCSF(4));
    qcsf.data.estAULCSF(qcsf.data.trial) = trapz(log10(qcsf.stimuli.frequency),qcsf.data.estSensitivity(qcsf.data.trial,:));
    
    if ~isfield(qcsf.data,'cResponses');
        qcsf.data.cResponses=[];
        qcsf.data.incResponses=[];
    end
    
    plotJitter = .02;
    plotData = [log10(qcsf.stimuli.nextFrequency)+plotJitter*randn -log10(qcsf.stimuli.nextContrast)+plotJitter*randn];
    
    if response,
        qcsf.data.cResponses=[plotData;
                              qcsf.data.cResponses];
    else
        qcsf.data.incResponses=[plotData;
                                qcsf.data.incResponses];
    end

function [Pc_output,H_output] = findCSF_Probability(LOGFREQ,LOGX,csfMode,csfParam,guessRate,lapseRate)

if strcmp(csfMode,'grid')
    qcsf = csfParam;
    
    [logGain,logCenter,octaveWidth,logTrunc] = ndgrid(qcsf.parameters.gain,qcsf.parameters.center,...
                                            qcsf.parameters.width,qcsf.parameters.trunc);
    
elseif strcmp(csfMode,'column')
    sCSF = csfParam;
    
    logGain = sCSF(:,1);
    logCenter =sCSF(:,2);
    octaveWidth = sCSF(:,3);
    logTrunc = sCSF(:,4);
    
elseif strcmp(csfMode,'point')
    QCSF=csfParam;
    
    logGain = QCSF(1);
    logCenter =QCSF(2);
    octaveWidth = QCSF(3);
    logTrunc = QCSF(4);
end
    
logTau = -findQCSF(LOGFREQ,logGain,logCenter,octaveWidth,logTrunc);

Pc_output = min(1-lapseRate,guessRate + (1-guessRate).*(1- exp(-10.^(2.*(LOGX-logTau)))));                 

H_output = -Pc_output.*log(Pc_output)-(1-Pc_output).*log(1-Pc_output);         

    
function estCSF = analyzePosterior(qcsf) 
%         
     posterior=qcsf.data.prior;
     
     g_hat = sum(qcsf.parameters.gain*sum(sum(posterior,3),4)); 
     c_hat = sum(qcsf.parameters.center*(sum(sum(posterior,3),4))');
     w_hat = sum(qcsf.parameters.width*squeeze(sum(sum(posterior,1),2)));     
     t_hat = sum(qcsf.parameters.trunc*squeeze(sum(sum(posterior,1),2))');   

     estCSF=[g_hat c_hat w_hat t_hat];   
      
function H = computeEntropy(P)

H = -P.*log(P)-(1-P).*log(1-P);


function qcsf=plotSimulation(qcsf,varargin)
    matlabVersion=version;
    if ~isfield(qcsf.simulation,'figure');
        
        if eval(matlabVersion(1))<7,
            figH=open('quickCSFv6.fig');
        else
            figH=open('quickCSF.fig');
        end
        
        figAxes=get(figH,'Children');
            tmpAxes=[8 1 7 4 6 5 2 3];
            
            for n=1:length(figAxes),
                qAxes(n)=figAxes(tmpAxes(n));
            end
            
            qcsf.simulation.figure.handles=qAxes;
    else
            qAxes = qcsf.simulation.figure.handles;
    end
    
    subplot(qAxes(1)); 
        set(gca,'Fontname','Myriad','FontSize',20)
        pause(.1)

        l(2)=plot(log10(qcsf.stimuli.frequency),qcsf.simulation.trueSensitivity,'k');hold on
        set(l(2),'lineWidth',4)


         if ~isempty(qcsf.data.cResponses) 
            p_tmp1=plot(qcsf.data.cResponses(:,1),qcsf.data.cResponses(:,2),'bo');
            set(p_tmp1,'MarkerSize',10,'Markerfacecolor',[0 0 1])
         end
         hold on
         if ~isempty(qcsf.data.incResponses)
            p_tmp2=plot(qcsf.data.incResponses(:,1),qcsf.data.incResponses(:,2),'ro');
            set(p_tmp2,'MarkerSize',10,'Markerfacecolor',[1 0 0])
         end

        l(1)=plot(log10(qcsf.stimuli.frequency),qcsf.data.estSensitivity(qcsf.data.trial,:),'go-'); 
        set(l(1),'lineWidth',4,'MarkerSize',8,'MarkerFaceColor',[0 1 0])

        axis([log10(qcsf.stimuli.frequency([1 end]))' -log10(qcsf.stimuli.contrast([end 1]))'])
        set(gca,'Yaxislocation','left','Xtick',log10([.5 1 2 5 10 20]),'XtickLabel',([.5 1 2 5 10 20]),'Ytick',log10([2 10 50 200 500 1000]),'Yticklabel',([2 10 50 200 500 1000]))
        xlabel('spatial frequency')
        ylabel('sensitivity')
        title('qCSF Simulation')
               set(gca,'Fontname','Myriad','FontSize',20)

        tH = text(log10(20),log10(600),['t=',num2str(qcsf.data.trial)]);
        set(tH,'Fontname','Myriad','FontSize',16)


        p(1)=plot(log10(50),0,'bo','MarkerfaceColor',[0 0 1]);
        p(2)=plot(log10(50),0,'ro','MarkerfaceColor',[1 0 0]);

        leg1=legend(p,'Correct','Incorrect',2);
        set(leg1,'Fontsize',16);

        hold off

        subplot(qAxes(2));

        set(gca,'Fontname','Myriad','FontSize',12)
        [contourX,contourY]=ndgrid(log10(qcsf.stimuli.frequency),log10(qcsf.stimuli.contrast));

        tmpEntropy = reshape(qcsf.stimuli.nextEntropy,size(contourX'))'./max(qcsf.stimuli.nextEntropy(:));

        pcolor(contourX,-contourY,tmpEntropy);
        shading interp;
        axis([log10(qcsf.stimuli.frequency([1 end]))' -log10(qcsf.stimuli.contrast([end 1]))'])
        hold on

        l(1) = plot(log10(qcsf.stimuli.frequency),qcsf.data.estSensitivity(qcsf.data.trial,:),'w-'); 

        set(l(1),'lineWidth',3)
        title('information gain ')
        axis([log10(qcsf.stimuli.frequency([1 end]))' -log10(qcsf.stimuli.contrast([end 1]))'])
        axis off
        hold off    

    subplot(qAxes(3));
        hold off;

        pcolor(qcsf.parameters.gain,qcsf.parameters.center,log10(marginalize(qcsf.data.prior,[3 4]))')
        shading interp


            l(1) = line([qcsf.simulation.trueCSF(1) qcsf.simulation.trueCSF(1)],([qcsf.parameters.center(1) qcsf.parameters.center(end)])); 
            l(2) = line([qcsf.parameters.gain(1) qcsf.parameters.gain(end)],[qcsf.simulation.trueCSF(2) qcsf.simulation.trueCSF(2)]); 
            set(l(1:2),'Color','w') 


        set(gca,'Ytick',log10([.2 .5 1 2 5 10 20]),'YtickLabel','',...
            'Xtick',log10([2 10 50 200 500]),'Xticklabel','','Yaxislocation','left','Fontsize',16,'Fontname','Myriad');

        ylabel('peak frequency');
        xlabel('peak gain')

        hold on
        priorSamples=1000;
         S = discretesamplev6(qcsf.data.prior,priorSamples);
         [sGain,sCenter,sWidth,sTrunc] = ind2sub(size(qcsf.data.prior),S);

        sCSF = [qcsf.parameters.gain(sGain)' qcsf.parameters.center(sCenter)' qcsf.parameters.width(sWidth)' qcsf.parameters.trunc(sTrunc)'] ;
            
        p1=plot(sCSF(:,1)+.1*randn(priorSamples,1),sCSF(:,2)+.1*randn(priorSamples,1),'w.');
        set(p1,'MarkerSize',3);
%             set(gca,'Fontname','Myriad','FontSize',12)

    subplot(qAxes(4));

        hold off
        pcolor((qcsf.parameters.width),qcsf.parameters.trunc,log10(marginalize(qcsf.data.prior,[1 2]))')
        shading interp


            l(3)=line(([qcsf.simulation.trueCSF(3) qcsf.simulation.trueCSF(3)]),[qcsf.parameters.trunc(1) qcsf.parameters.trunc(end)]); 
            l(4) = line(([qcsf.parameters.width(1) qcsf.parameters.width(end)]),[qcsf.simulation.trueCSF(4) qcsf.simulation.trueCSF(4)]); 
            set(l(3:4),'Color','w') 


        hold on
        p2=plot(sCSF(:,3)+.05*randn(priorSamples,1),sCSF(:,4)+.05*randn(priorSamples,1),'w.');
        set(p2,'MarkerSize',3);
        set(gca,'Fontname','Myriad','FontSize',16,'Yticklabel','','Xticklabel','','xaxislocation','top')

        ylabel('truncation','rotation',90)
        xlabel('bandwidth')

            
        %%%%%%%%%%%%%%%%%%%%%%  PLotting marginals
        set(qAxes(5:8),'Fontname','Myriad','FontSize',12)

        subplot(qAxes(5));
            plot(qcsf.parameters.gain,marginalize(qcsf.data.prior,[2 3 4])./diff(qcsf.parameters.gain(1:2)));
            set(gca,'Xlim',qcsf.parameters.gain([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([2 10 50 200 500]),'Xticklabel',([2 10 50 200 500]),'xaxislocation','top')
            ylabel('pdf')
            
        subplot(qAxes(6));
            plot(marginalize(qcsf.data.prior,[1 3 4])./diff(qcsf.parameters.center(1:2)),qcsf.parameters.center);
            set(gca,'Ylim',qcsf.parameters.center([1 end]),'Xtick',[0 5],'Xlim',[0 8],'Ytick',log10([.2 .5 1 2 5 10 20]),'YtickLabel',[.2 .5 1 2 5 10 20],'yaxislocation','right','xaxislocation','top')
            ylabel('cycles per degree ')
            xlabel('pdf')
        
        subplot(qAxes(8));
            plot(qcsf.parameters.width,marginalize(qcsf.data.prior,[1 2 4])./diff(qcsf.parameters.trunc(1:2)));
            set(gca,'Xlim',qcsf.parameters.width([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([1 2 4 8]),'XtickLabel',([1 2 4 8]))
            ylabel('pdf')
            xlabel('octaves ')
            
        subplot(qAxes(7));
            plot(marginalize(qcsf.data.prior,[1 2 3])./diff(qcsf.parameters.trunc(1:2)),qcsf.parameters.trunc);
            set(gca,'Ylim',qcsf.parameters.trunc([1 end]),'Xlim',[0 8],'Xtick',[0 5],'xaxislocation','bottom','yaxislocation','right','ytick',log10([.02 .05 .2 .5 2]),'ytickLabel',[.02 .05 .2 .5 2])
            xlabel('pdf')
            ylabel('decimal log units ')
   
            
        qcsf.stimuli=rmfield(qcsf.stimuli,{'nextContrast' 'nextFrequency' 'nextEntropy'});
         
%                
function qcsf=plotExperiment(qcsf)
       
    matlabVersion=version;
        
    if eval(matlabVersion(1))<7,
        figH=open('quickCSFv6.fig');
    else
        figH=open('quickCSF.fig');
    end
    
    figAxes=get(figH,'Children');
    tmpAxes=[8 1 7 4 6 5 2 3];

    for n=1:length(figAxes),
        qAxes(n)=figAxes(tmpAxes(n));
    end
     
    subplot(qAxes(1)); 
        set(gca,'Fontname','Myriad','FontSize',20)
        pause(.1)

         if ~isempty(qcsf.data.cResponses) 
            p_tmp1=plot(qcsf.data.cResponses(:,1),qcsf.data.cResponses(:,2),'bo');
            set(p_tmp1,'MarkerSize',10,'Markerfacecolor',[0 0 1])
         end
         hold on
         if ~isempty(qcsf.data.incResponses)
            p_tmp2=plot(qcsf.data.incResponses(:,1),qcsf.data.incResponses(:,2),'ro');
            set(p_tmp2,'MarkerSize',10,'Markerfacecolor',[1 0 0])
         end

        l(1)=plot(log10(qcsf.stimuli.frequency),qcsf.data.estSensitivity(qcsf.data.trial,:),'go-'); 
        set(l(1),'lineWidth',4,'MarkerSize',8,'MarkerFaceColor',[0 1 0])

        axis([log10(qcsf.stimuli.frequency([1 end]))' -log10(qcsf.stimuli.contrast([end 1]))'])
        set(gca,'Yaxislocation','left','Xtick',log10([.5 1 2 5 10 20]),'XtickLabel',([.5 1 2 5 10 20]),'Ytick',log10([2 10 50 200 500 1000]),'Yticklabel',([2 10 50 200 500 1000]))
        xlabel('spatial frequency')
        ylabel('sensitivity')
        title('qCSF Summary')
               set(gca,'Fontname','Myriad','FontSize',20)

        tH = text(log10(20),log10(600),['t=',num2str(qcsf.data.trial)]);
        set(tH,'Fontname','Myriad','FontSize',16)

        p(1)=plot(log10(50),0,'bo','MarkerfaceColor',[0 0 1]);
        p(2)=plot(log10(50),0,'ro','MarkerfaceColor',[1 0 0]);

        leg1=legend(p,'Correct','Incorrect');
        set(leg1,'Fontsize',16);

        hold off

        subplot(qAxes(3));
        hold off;

        pcolor(qcsf.parameters.gain,qcsf.parameters.center,log10(marginalize(qcsf.data.prior,[3 4]))')
        shading interp

        set(gca,'Ytick',log10([.2 .5 1 2 5 10 20]),'YtickLabel','',...
            'Xtick',log10([2 10 50 200 500]),'Xticklabel','','Yaxislocation','left','Fontsize',16,'Fontname','Myriad');

        ylabel('peak frequency');
        xlabel('peak gain')

        hold on
        priorSamples=1000;
         S = discretesamplev6(qcsf.data.prior,priorSamples);
         [sGain,sCenter,sWidth,sTrunc] = ind2sub(size(qcsf.data.prior),S);

        sCSF = [qcsf.parameters.gain(sGain)' qcsf.parameters.center(sCenter)' qcsf.parameters.width(sWidth)' qcsf.parameters.trunc(sTrunc)'] ;
            
        p1=plot(sCSF(:,1)+.1*randn(priorSamples,1),sCSF(:,2)+.1*randn(priorSamples,1),'w.');
        set(p1,'MarkerSize',3);
%             set(gca,'Fontname','Myriad','FontSize',12)

    subplot(qAxes(4));

        hold off
        pcolor((qcsf.parameters.width),qcsf.parameters.trunc,log10(marginalize(qcsf.data.prior,[1 2]))')
        shading interp

        hold on
        p2=plot(sCSF(:,3)+.05*randn(priorSamples,1),sCSF(:,4)+.05*randn(priorSamples,1),'w.');
        set(p2,'MarkerSize',3);
        set(gca,'Fontname','Myriad','FontSize',16,'Yticklabel','','Xticklabel','','xaxislocation','top')

        ylabel('truncation','rotation',90)
        xlabel('bandwidth');
        
        for n=1:size(qcsf.data.estCSF,1),
            estAULCSF(n) = trapz(log10(qcsf.stimuli.frequency),qcsf.data.estSensitivity(n,:));
        end
        
        subplot(qAxes(2))
        plot(log(1:length(estAULCSF)),log(estAULCSF),'k')

        set(gca,'Xaxislocation','top','Xtick',log([3 10 30 100 300]),'Xticklabel',([3 10 30 100 300]),'Yaxislocation','right','FontSize',8,...
         'Ytick',log([1 2 4 8]),'Yticklabel',[1 2 4 8]);

        ylabel('AULCSF')
        axis([0 log(qcsf.data.trial) log(1) log(10)]);
        text(log(1),log(8),['AULCSF = ' sprintf('%4.4f', estAULCSF(qcsf.data.trial))]) 
        
    %%%%%%%%%%%%%%%%%%%%%%  PLotting marginals of the posteriors
        set(qAxes(5:8),'Fontname','Myriad','FontSize',12)

        subplot(qAxes(5));
            plot(qcsf.parameters.gain,marginalize(qcsf.data.prior,[2 3 4])./diff(qcsf.parameters.gain(1:2)));
            set(gca,'Xlim',qcsf.parameters.gain([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([2 10 50 200 500]),'Xticklabel',([2 10 50 200 500]),'xaxislocation','top')
            ylabel('pdf')
            
        subplot(qAxes(6));
            plot(marginalize(qcsf.data.prior,[1 3 4])./diff(qcsf.parameters.center(1:2)),qcsf.parameters.center);
            set(gca,'Ylim',qcsf.parameters.center([1 end]),'Xtick',[0 5],'Xlim',[0 8],'Ytick',log10([.2 .5 1 2 5 10 20]),'YtickLabel',[.2 .5 1 2 5 10 20],'yaxislocation','right','xaxislocation','top')
            ylabel('cycles per degree ')
            xlabel('pdf')
        
        subplot(qAxes(8));
            plot(qcsf.parameters.width,marginalize(qcsf.data.prior,[1 2 4])./diff(qcsf.parameters.trunc(1:2)));
            set(gca,'Xlim',qcsf.parameters.width([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([1 2 4 8]),'XtickLabel',([1 2 4 8]))
            ylabel('pdf')
            xlabel('octaves ')
            
        subplot(qAxes(7));
            plot(marginalize(qcsf.data.prior,[1 2 3])./diff(qcsf.parameters.trunc(1:2)),qcsf.parameters.trunc);
            set(gca,'Ylim',qcsf.parameters.trunc([1 end]),'Xlim',[0 8],'Xtick',[0 5],'xaxislocation','bottom','yaxislocation','right','ytick',log10([.02 .05 .2 .5 2]),'ytickLabel',[.02 .05 .2 .5 2])
            xlabel('pdf') 
            ylabel('decimal log units ')
            
%%%%%%%%%%%%%%%%%%%%%%  Plot the experiment's initial marginal priors for comparison
        
        subplot(qAxes(5));
        hold on
            plot(qcsf.parameters.gain,marginalize(qcsf.data.prior0,[2 3 4])./diff(qcsf.parameters.gain(1:2)),'r-');
            set(gca,'Xlim',qcsf.parameters.gain([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([2 10 50 200 500]),'Xticklabel',([2 10 50 200 500]),'xaxislocation','top')
            ylabel('pdf')
            
        subplot(qAxes(6));
        hold on
            plot(marginalize(qcsf.data.prior0,[1 3 4])./diff(qcsf.parameters.center(1:2)),qcsf.parameters.center,'r-');
            set(gca,'Ylim',qcsf.parameters.center([1 end]),'Xtick',[0 5],'Xlim',[0 8],'Ytick',log10([.2 .5 1 2 5 10 20]),'YtickLabel',[.2 .5 1 2 5 10 20],'yaxislocation','right')
            ylabel('cycles per degree ')
            xlabel('pdf')
        
        subplot(qAxes(8));
        hold on;
            plot(qcsf.parameters.width,marginalize(qcsf.data.prior0,[1 2 4])./diff(qcsf.parameters.trunc(1:2)),'r-');
            set(gca,'Xlim',qcsf.parameters.width([1 end]),'Ytick',[0 5],'Ylim',[0 8],'Xtick',log10([1 2 4 8]),'XtickLabel',([1 2 4 8]))
            ylabel('pdf')
            xlabel('octaves ')
            
        subplot(qAxes(7));
        hold on;
            plot(marginalize(qcsf.data.prior0,[1 2 3])./diff(qcsf.parameters.trunc(1:2)),qcsf.parameters.trunc,'r-');
            set(gca,'Ylim',qcsf.parameters.trunc([1 end]),'Xlim',[0 8],'Xtick',[0 5],'xaxislocation','top','yaxislocation','right','ytick',log10([.02 .05 .2 .5 2]),'ytickLabel',[.02 .05 .2 .5 2])
            xlabel('pdf')
            ylabel('decimal log units ')
     
       


%%

       


