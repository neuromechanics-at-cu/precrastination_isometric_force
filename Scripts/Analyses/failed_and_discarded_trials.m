
%% Load Data
if ~exist('M')
    load('../../Data/Processed Data/processed_data.mat','M')
end
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end


%% Failed Trials 
failed_choice_stats = [];
failed_fam_stats = [];
for ii = 1:length(subjs)
    for kk = 1:length(M{ii}.failedTrials)
        failed_fam_stats = [failed_fam_stats;...
                        ii,M{ii}.failedTrials(kk)]
    end
    for kk = 1:length(Mchoice{ii}.failedTrials)
        failed_choice_stats = [failed_choice_stats;...
                               ii,Mchoice{ii}.failedTrials(kk)]
    end
end
% No failed trials

%% Thrown Out Trials
% These are due to issues with the robot / code.
% Issues arise consistently after rest periods, 
% The other is less clear, but some choice periods are very short
% (<0.2sec), trials here are inconsistent.
badtrial_info_A = [];% Hills
badtrial_info_B = [];% Valleys
for ii = 1:length(subjs)
    for jj = 1:910
        if isnan(Mchoice{ii}.choice(jj))
            if Mchoice{ii}.trialnumber(jj) <= 455 %is hill
                badtrial_info_A = [badtrial_info_A;
                                 ii,jj];
            else
                badtrial_info_B = [badtrial_info_B;
                                 ii,jj];
            end
        end    
    end
end