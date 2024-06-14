% Get fixation data
function [eyeData] = ImportFixationData(displayFps, DATA_DIR, save_flag, jitter_flag)

if nargin <=2 
    save_flag = 0;
    jitter_flag = 0;
end
% DATA_DIR        = '../Data/eyemovements/';
if ~exist(DATA_DIR)
    DATA_DIR        = '/home/vaishnavi/UW/Sight_Restoration/p2p_optoExpt/Data/eyemovements/';
    if ~exist(DATA_DIR)
        sprintf("Please enter a valid directory containing eye-movement data. Thank you!")
        return
    end
end

files = dir(fullfile(DATA_DIR, '*.mat'));
aggr_fix_dur = [];

% Display parameters; First row corresponds to original values as given in
% DOVES eye movement database
% Ref: http://live.ece.utexas.edu/publications/2009/ivl_sv_feb09.pdf
% second row corresponds to our display's parameters
% All dimensions are in cm and angles in degrees



org.fps = 200;
org.res = [1024, 768];
org.depth = 134;
org.span = [17, 13];
org.dims = [];
org.pixperdeg = round(org.res(1)/org.span(1));

disp.fps = displayFps; % NOTE: change our fps from 60 to 120 when using the lab monitor instead of my dev-machine
disp.res = [1920 1080];
disp.depth = 80;
disp.dims = [71, 40];
disp.span = PixelToVisAngle([disp.dims, disp.depth], 0);
disp.pixperdeg(:) = round(disp.res(1)./disp.span(1));

scaleFac = [(disp.pixperdeg/org.pixperdeg)*[1,1], disp.fps/org.fps];

%Load the fixations from all subjects for all images.
for img_id = 1:length(files) % iterating through images
    load (append(DATA_DIR, files(img_id).name)); %loads subj_names_list, fix_data, eye_data
    for sub_id = 1:length(fix_data) % iterating through subjects
        if jitter_flag ~= 0
            eyeData{1, sub_id} = [];
        end
        %         if img_id == 1
        %             fix_dur.("sub_" + string(sub_id)) = [];
        %             fix_loc.("sub_" + string(sub_id)) = [];
        %         end
        %         temp_dur = fix_data{sub_id}(3,:)/fps;
        %         temp_dur(temp_dur < 0) = 0;
        %         temp_dur = temp_dur(~isnan(temp_dur));
        %         [temp_theta, temp_rho] = cart2pol(fix_data{sub_id}(1,:), fix_data{sub_id}(2,:));
        %
        %         fix_dur.("sub_" + string(sub_id)) = [fix_dur.("sub_" + string(sub_id)), temp_dur] ;
        %         fix_loc.("sub_" + string(sub_id)) = [fix_loc.("sub_" + string(sub_id)), [temp_theta; temp_rho]];
        %         aggr_fix_dur = [aggr_fix_dur, temp_dur];

        temp = [];
        temp = diff(fix_data{1, sub_id},1,2);
        if size(temp,2) > 1
            temp(1:2,1) = 0;
%             temp(1,:) = ((temp(1,:)*scaleFac(1)));
%             temp(2,:) = (temp(2,:)*scaleFac(2));
            temp(1,2:end) = ((temp(1,2:end)*scaleFac(1)));
            temp(2,2:end) = (temp(2,2:end)*scaleFac(2));
        end

        if jitter_flag == 1
            temp(3,:) = floor((fix_data{1, sub_id}(3,1:end-1)*scaleFac(3))/10);
            temp(4,:) = FindProjection(temp(1:2,:), -1);
            temp(5,:) = FindProjection(temp(1:2,:), +1);
            temp = repmat(temp(:,2:end), 1, 30);
            temp(1:2,1) = 0; temp(4:5,1) = 0;
            eyeData{1, sub_id} = [eyeData{1, sub_id}, temp];
        elseif(jitter_flag) == 2
            eyeData{1, sub_id} = [0;0;300;0;0] ;
            
        else
            temp(3,:) = floor((fix_data{1, sub_id}(3,1:end-1)*scaleFac(3))/2);
            temp(4,:) = FindProjection(temp(1:2,:), 1);
            temp(5,:) = FindProjection(temp(1:2,:),-1);
            eyeData{img_id, sub_id} = temp;
        end

        %         figure(1)
        %         plot(fix_loc.("sub_" + string(sub_id))(1,:), fix_loc.("sub_" + string(sub_id))(2,:),'ko', 'MarkerFaceColor','g');
        %         hold on
        %         plot(fix_loc.("sub_" + string(sub_id))(1,:), fix_loc.("sub_" + string(sub_id))(2,:));
        %         hold off
        %         axis ij;
    end
end

% Fitting a poisson distribution for each subject to check if there is any variation
% across subjects
% no_subs = length(fieldnames(fix_dur));
% param = zeros(no_subs,1);
% for sub_id = 1:no_subs
%         param(sub_id) = poissfit(fix_dur.("sub_" + string(sub_id)));
% end
%
% % Doesn't seem like a significant difference
% est_param_clt = mean(param);
% est_param_aggr = poissfit(aggr_fix_dur);
%
% figure(2)
% subplot(1,2,1)
% histfit(aggr_fix_dur, 100, 'poisson')
% xlabel("Fixation durations across all subjects")
% subplot(1,2,2)
% histfit(param)
% xlabel("Lambda distribution for a poisson fit across individual subjects")

% Fitting a 2-D gaussian distribution for the fixation locations
% for sub_id = 1:no_subs
%     figure(3)
%     subplot(1,2,1)
%     histfit(fix_loc.("sub_" + string(sub_id))(2,:))
%     xlabel("Amplitude of change in fixation location")
%     subplot(1,2,2)
%     histfit(fix_loc.("sub_" + string(sub_id))(1,:))
%     xlabel("Angle of change in fixation location")
% end

    if save_flag
        save("eye_data.mat", "eyeData");
    end
end
%%