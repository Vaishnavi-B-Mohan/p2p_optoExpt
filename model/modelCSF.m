classdef modelCSF

    properties (Access = public)
        subjectID;
        modelData = [];
        % initial (default) D.H. Kelly parameters
        initParams = struct('param1', 6.1, 'param2', 7.3, 'param3', 45.9, 'b', 2);
    end

    methods (Access = public)
        function this = modelCSF(subjectID)
            this.subjectID = subjectID;
        end

        function this = collateData(this, subjectID, save_flag)
            fileList = dir(['../' filesep 'output' filesep subjectID filesep subjectID '_' '*.mat']);
            % ignore all files except experimental conditions
            indices= cellfun(@(x) contains(x,'DEMO') || ...
                contains(x,'error'),{fileList.name}, 'Unif',0);
            fileList = fileList(~[indices{:}]);

            %             dataTable = table;
            %             varTypes= {'string', 'string','double','double','double','double'};
            %             varNames = {'SubID', 'Condition','TF','SF', 'Contrast', 'Response'};

            baseline = []; opto = []; eye = [];
            for i= 1:length(fileList)
                %% load a single file
                fname = fileList(i).name;
                splitfname = strsplit(fname, '_');
                condition = splitfname{6};
                file(i) = load([fileList(i).folder filesep fname]);
                data = file(i).data; keys = data.keys;

                for idx = 1: data.numEntries
                    temp = strsplit(keys(idx), '_');
                    FF = str2double(temp(2));
                    FF_data = data(keys(idx));
                    numTrials = FF_data.qcsf.data.trial;
                    iter_dataMat = zeros(numTrials, 4);
                    %                     set(gcf,'Renderer','painters');
                    if strcmp(temp(1), 'TF')
                        iter_dataMat(:,1) = FF;
                        iter_dataMat(:,2) = FF_data.qcsf.data.history(:,2);
                        iter_dataMat(:,3) = log10(ones(numTrials, 1) ./ FF_data.qcsf.data.history(:,3)); % Sensitivity = log10(1/Contrast)
                        iter_dataMat(:,4) = FF_data.qcsf.data.history(:,4);

                        %                         correctIdx = (iter_dataMat(:,4) == 1);
                        %                         scatter(iter_dataMat(correctIdx,2), iter_dataMat(correctIdx, 3), 'filled');

                    elseif strcmp(temp(1), 'SF')
                        iter_dataMat(:,1) = FF_data.qcsf.data.history(:,2);
                        iter_dataMat(:,2) = FF;
                        iter_dataMat(:,3) = log10(ones(numTrials, 1) ./ FF_data.qcsf.data.history(:,3)); % Sensitivity = log10(1/Contrast)
                        iter_dataMat(:,4) = FF_data.qcsf.data.history(:,4);

                        %                         correctIdx = (iter_dataMat(:,4) == 1);
                        %                         scatter(iter_dataMat(correctIdx,1), iter_dataMat(correctIdx, 3), 'filled');
                    end
                    %                     title(keys(idx), 'Interpreter','none')
                    %                     colorbar;
                    %                     grid on
                    %                     ylabel('Sensitivity')
                    if strcmp(condition, 'baseline')
                        baseline = [baseline; iter_dataMat]; %% col1 == TF; col2 == SF; col3 == contrast; col4 == response
                    elseif strcmp(condition, 'opto')
                        opto = [opto; iter_dataMat];
                    elseif strcmp(condition, 'eye')
                        eye = [eye; iter_dataMat];
                    end
                end

            end
              
            if ~isempty(baseline)
                exptData(1).condition = 'baseline';
                exptData(1).data = baseline; 
            end
            if ~isempty(opto)
                exptData(2).condition = 'opto';
                exptData(2).data = opto; 
            end
            if ~isempty(eye)
                exptData(3).condition = 'eye';
                exptData(3).data = eye;
            end
            

            if save_flag
                collatedData_File = string(['../analysis_data' filesep subjectID '_collatedData.mat']);
                save(collatedData_File, "exptData");
%                 baseline_file = string(['../analysis_data' filesep subjectID '_baseline.mat']);
%                 save(baseline_file, "baseline");
%                 opto_file = string(['../analysis_data' filesep subjectID '_opto.mat']);
%                 save(opto_file, "opto");
%                 eye_file = string(['../analysis_data' filesep subjectID '_eye.mat']);
%                 save(eye_file, "eye");
            else
                
                this.modelData(1).condition = 'baseline'; this.modelData(2).condition = 'opto'; this.modelData(3).condition = 'eye';
                this.modelData(1).data = baseline; this.modelData(2).data = opto; this.modelData(3).data = eye;
                plotcolors = ['b'; 'r'; 'g'];
                markerfacecolors = [0, 0, 1; 1, 0, 1; 1, 1, 0];
                for i = 1:length(this.modelData)
%                     subplot(length(this.modelData), 1, i)
                    nTrials = size(this.modelData(i).data,1);
                    this.modelData(i).data = [this.modelData(i).data, zeros(nTrials,1)];
                    [this.modelData(i).data(:,end), this.modelData(i).bestParams] = fitMLE(this, this.modelData(i).data);
                    [plotsf,plottf] = meshgrid(2.^linspace(0.5,5,101));
                    %                     surf(log(plotsf(1,:)),log(plottf(:,1)), this.modelData(i).bestParams.sens , 'FaceColor', plotcolors(i), 'FaceAlpha',.2,'EdgeColor','none');
                    surf(log10(plotsf(1,:)),log10(plottf(:,1)), this.modelData(i).bestParams.sens , 'FaceColor', plotcolors(i), 'FaceAlpha',.3,'EdgeColor','none');
%                     logx2raw(10)
%                     logy2raw(10)
                    %                     surf(log(plotsf(1,:)),log(plottf(:,1)), initParams.sens , 'FaceColor', plotcolors(i), 'FaceAlpha',.2,'EdgeColor','none');
                    set(gca,'Zlim',[0,2.5]);
                    xlabel('SF');
                    ylabel('TF');
                    grid on;
                    view(25,25)
                    zlabel('Sens')
%                     hold on;
%                     scatter3(log10(this.modelData(i).data(:,2)), log10(this.modelData(i).data(:,1)), this.modelData(i).data(:,5), 'MarkerEdgeColor','k','MarkerFaceColor',markerfacecolors(i,:));
                    logx2raw(10)
                    logy2raw(10)
                                        legend('baseline', 'opto', 'eye');
%                     legend('fit', 'data');
                    hold on;
                end
               

                title(subjectID,'Interpreter','none');
            end
        end


        function [Sens, bestParams] = fitMLE(this, inputData)
            [plotsf,plottf] = meshgrid(2.^linspace(0.5,5,101));
            this.initParams.sens = calculateKelly(this, this.initParams,plotsf,plottf);

            tf =  inputData(:,1);
            sf =  inputData(:,2);
            contrast =  0.1.^inputData(:,3);
            resp = inputData(:,4);

            % Find the sensitivity for each stimulus
            S = calculateKelly(this, this.initParams,sf,tf);

            % convert from log10 sensitivity to threshold contrast.  Save them as
            % threshold values for the Weibull function
            this.initParams.t = .1.^S;    % or 1./S.^10
            this.initParams.err = fitKelly(this, this.initParams,sf(:),tf(:),contrast(:),resp);

            freeParams = {'param1','param2','param3'};
            bestParams = fit('fitKelly', this.initParams,freeParams,sf(:),tf(:),contrast(:),resp);
            bestParams.err = fitKelly(this, bestParams,sf(:),tf(:),contrast(:),resp);
            bestParams.sens = calculateKelly(this, bestParams,plotsf,plottf);
            Sens = -log(inputData(:,3));
        end


        function [err,prob] = fitKelly(this, p,sf,tf,c,resp)
            % [err,prob] = fitKelly(p,sf,tf,c,resp)
            %
            % calculates the -log likelihood of data in 'resp' of the TCSF model from
            % DH Kelly (1979) and the weibull psychometric fuction.

            % Calculate the TCSF surface
            S = calculateKelly(this,p,sf,tf);

            % Convert from log10 sensitivity to contrast thresholds
            p.t = .1.^S(:);
            %             p.t = 10.^(-S(:));

            % Calculate P(correct) for each trial
            prob = weibull(this,p,c(:));
            prob = min(prob,.99); % standard hack to avoid log(0) NaNs

            % -log likelihood calculation
            err = -sum(resp.*log(prob) + (1-resp).*log(1-prob));

        end

        function sens = calculateKelly(this, p, sf, tf)
            % sens = calculateKelly(p, sf, tf)
            %
            % 4/5/2024  gmb  removed adding NaN's to low senstivity values

            % calculate Kelly values
            alpha = 2*pi*sf;	%used in the formula, related to spatial frequency
            v = tf ./ sf;
            k = p.param1 + p.param2 *(abs(log10( v / 3 ))) .^ 3; %p.param1 = 6.1; p.param2 = 7.3;
            amax = p.param3 ./ (v + 2); % p.param3 = 45.9;
            sens = k .* v .* (alpha .* alpha) .* exp(-2*alpha ./ amax);

            %	Pull out the ones above 2;
            sens(sens<2) = 2;% ones(size(sens(l)))*(0/0);  %	This sets the bad ones to NaN
            sens = log10(sens/2);
        end


        function [p] = weibull(this,params, x)
            % [p] = Weibull(params, x)
            %
            % The Weibull function based on this equation:
            %
            % k = (-log((1-e)/(1-g)))^(1/b)
            % f(x) = 1 - ((1-g) * exp(-(k*x/t).^b))
            %
            % Where g is performance expected at chance, e is performance level that
            % defines the threshold, b is the slope of the Weibull function, and t is
            % the threshold
            %
            % Inputs:
            %   params      A structure containing the parameters of the Weibull
            %               function:
            %       b       Slope
            %       t       Stimulus intensity threshold as defined by 'params.e'.
            %               When x = 'params.t', then y = 'params.e'
            %       g       Performance expected at chance, proportion
            %       e       Threshold performance, proportion
            %
            %   x           Intensity values of the stimuli
            %
            % Output:
            %   p           Output of the Weibull function as a function of the
            %               intensity values, x

            % Written by G.M. Boynton - 11/13/2007
            % Edited by Kelly Chang - February 13, 2017
            % Edited by Ione Fine - February 22, 2017

            %% Evaluate Weibull Function
            if ~isfield(params, 'g')
                params.g = 0.8;
            end
            if ~isfield(params, 'e')
                params.e = (0.5)^(1/3);
            end

            k = (-log((1-params.e)/(1-params.g)))^(1/params.b);
            p = 1 - ((1-params.g).* exp(-(k*x./params.t).^params.b));
        end


    end
end
