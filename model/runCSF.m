clear all; close all; clc;

subjectID = dir(['..' filesep 'output' filesep 'NS_*']);
subjectID(strncmp({subjectID.name}, '.', 1)) = []; %remove hidden files

figure(1)
for i = 1:size(subjectID,1)
%     figure(i)
    disp(subjectID(i).name)
    model = modelCSF(subjectID(i).name);
    subplot(3,3,i)
    model = model.collateData(subjectID(i).name, 0); % Second parameter takes save_flag. If set to 1, it saves the collated files

end

print('Lalala')