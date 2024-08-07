%% Construct a GLME for Choice Prefences
% TABLE 1 in the manuscript


%%%  Fixed Effects:
% Hill or Valley 
% Difference in Hill/Valley Duration 
% Difference in Hill/Valley Onset Time
%
%%% Random Effects
% Subject
%
%%% Model Type
% logit(choice) ~ effects

%% Load Data
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end

%% Make A Table
choice_table = [];
for jj = 1:length(subjs)
    for kk = 1:length(Mchoice{jj}.choice)
        if ~isnan(Mchoice{jj}.choice(kk))
%             choice_table = [choice_table;
%                             %subj   hill/vall                       diff duration                                               diff onset                                             chose left ...
% %                             jj,     Mchoice{jj}.trialnumber(kk) <= 455, Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk), -1*Mchoice{jj}.choice(kk)+2];
% %                             subj   hill/vall                       diff duration                                               diff onset                                             chose left ...
%                             jj,     Mchoice{jj}.trialnumber(kk) <= 455, Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk), -1*Mchoice{jj}.choice(kk)+2,];
% 
          choice_table = [choice_table;
%                             subj   hill/vall                             % First half/Second Half                 diff duration                                               diff onset                                             chose left ...
%                             jj,     Mchoice{jj}.trialnumber(kk) <= 455, Mchoice{jj}.trialnumber_bysession(kk)<227, Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk), -1*Mchoice{jj}.choice(kk)+2,];
                            jj,     Mchoice{jj}.trialnumber(kk) <= 455, Mchoice{jj}.trialnumber_bysession(kk), Mchoice{jj}.trialnumber_bysession(kk)<227, Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk), -1*Mchoice{jj}.choice(kk)+2,];
        end
    end
end

% Define Variable Types and Name Variables
% choice_table = array2table(choice_table,'VariableNames',{'subj','hill_or_valley','delta_duration','delta_onset','chose_left'});
% choice_table = array2table(choice_table,'VariableNames',{'subj','hill_or_valley','first_or_second_half','delta_duration','delta_onset','chose_left'});
choice_table = array2table(choice_table,'VariableNames',{'subj','hill_or_valley','trial_num','first_or_second_half','delta_duration','delta_onset','chose_left'});
choice_table.subj = categorical(choice_table.subj);
choice_table.hill_or_valley = categorical(choice_table.hill_or_valley);
choice_table.first_or_second_half = categorical(choice_table.first_or_second_half);
choice_table.chose_left = logical(choice_table.chose_left);

%% Fit Model
% mdl = fitglme(choice_table,'chose_left ~ hill_or_valley*delta_duration*delta_onset + (1|subj)','Distribution','Binomial','Link','Logit')
% mdl = fitglme(choice_table,'chose_left ~ hill_or_valley*delta_duration*delta_onset*first_or_second_half + (1|subj)','Distribution','Binomial','Link','Logit')

% mdl = fitglme(choice_table,'chose_left ~ hill_or_valley*delta_duration*delta_onset+first_or_second_half*hill_or_valley + (1|subj)','Distribution','Binomial','Link','Logit')

mdl = fitglme(choice_table,'chose_left ~ hill_or_valley*delta_duration*delta_onset+trial_num*hill_or_valley + (1|subj)','Distribution','Binomial','Link','Logit')