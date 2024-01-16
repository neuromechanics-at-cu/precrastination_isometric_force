%% Load in Data
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('M')
    load('../../Data/Processed Data/processed_data.mat','M')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end



%% Familiarization Phase Accuracy
% Compile Accuracy
acc_vec = [];
for jj = 1:length(subjs)
    acc_vec = [acc_vec,M{jj}.accuracy'];
end
% Calculate Mean and SE, in % (not px)
nanmean(acc_vec/824)
nanstd(acc_vec/824)/sqrt(length(subjs))


%% Choice Phase Accuracy
% Compile Accuracy
acc_vec = [];
for jj = 1:length(subjs)
    acc_vec = [acc_vec,Mchoice{jj}.accuracy'];
end
% Calculate Mean and SE, in % (not px)
nanmean(acc_vec/824)
nanstd(acc_vec/824)/sqrt(length(subjs))


%% Linear Models 

% Construct Table of Choices
acc_table = [];
for jj = 1:length(subjs)
    for kk = 1:910
        if ~isnan(Mchoice{jj}.choice(kk))
        % Subject // Delta Accuracy // choice
        acc_table = [acc_table;
                jj,(M{jj}.mean_accuracytot(Mchoice{jj}.Ltrialtype(kk)+1) - M{jj}.mean_accuracytot(Mchoice{jj}.Rtrialtype(kk)+1))/824,...
                Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk),-1*Mchoice{jj}.choice(kk)+2];
        end
    end
end
acc_table = array2table(acc_table,'VariableNames',{'subj','delta_accuracy','delta_duration','delta_onset','chose_left'});
acc_table.subj = categorical(acc_table.subj);
acc_table.chose_left = logical(acc_table.chose_left);


% Make new table of Just unique combos
nComb = nchoosek(1:14,2);
nComb =[nComb;nComb+14];
% Hack dOnset and dDuration Table
for ii = 1:size(nComb,1)
    ind = find(M{1}.trialtype == nComb(ii,1)-1,1);
    topdur(ii,1) = M{1}.topduration(ind);
    onsets(ii,1) = M{1}.toponset(ind);
    ind2 = find(M{1}.trialtype == nComb(ii,2)-1,1);
    topdur(ii,2) = M{1}.topduration(ind2);
    onsets(ii,2) = M{1}.toponset(ind2);
end

acc_table_combos = [];
for jj = 1:length(subjs)
    for kk = 1:size(nComb,1)
        % Subject // Delta Accuracy // choice
        acc_table_combos = [acc_table_combos;
                jj,(M{jj}.mean_accuracytot(nComb(kk,1)) - M{jj}.mean_accuracytot(nComb(kk,2)))/824,...
                topdur(kk,1) - topdur(kk,2),onsets(kk,1) - onsets(kk,2),...
                kk<=91];
    end
end
acc_table_combos = array2table(acc_table_combos,'VariableNames',{'subj','delta_accuracy','delta_duration','delta_onset','session_type'});
acc_table_combos.subj = categorical(acc_table_combos.subj);
acc_table_combos.session_type = categorical(acc_table_combos.session_type);


%% Diff in Accuracy and Preference
% Choice Model
fitglme(acc_table,'chose_left ~ delta_accuracy + (1|subj)','Distribution','Binomial','Link','Logit')


%% Diff in Accuracy and Diff in Duration and Onset
fitlme(acc_table_combos,'delta_accuracy ~ delta_duration + delta_onset + (1|subj)')
