%% Plot Linear Model Results
% Figure 3 in the Manuscript


%% Load in Data
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end

%% Initialize Variables
session_strings = {'Hill','Valley'};
unique_delta_duration = -4:2:4; 
unique_delta_onset = -10:2:10;
unique_delta_offset = -10:2:10;
pChoice_delta_duration = nan(length(subjs),2,length(unique_delta_duration));
pChoice_delta_onset = nan(length(subjs),2,length(unique_delta_onset));
pChoice_delta_offset = nan(length(subjs),2,length(unique_delta_offset));

colorvec = [0,0,100/255;
         0,100/255,0];

%Which choices do you want to visualize? 
nChoice = 1;%or later 
% Example:
% 1: all 5 choices
% 3: last 3 choices

%% Construct Massive Table
TABLE = [];
for ii = 1:length(subjs)
        delta_duration  = Mchoice{ii}.Ltopduration - Mchoice{ii}.Rtopduration;
        delta_onset  = Mchoice{ii}.Ltoponset - Mchoice{ii}.Rtoponset;
        delta_offset = delta_duration + delta_onset;
        trialnumber = Mchoice{ii}.trialnumber;
        choicetypecount = Mchoice{ii}.choicetypecount;
        choseLeft = -1*Mchoice{ii}.choice+2;
        subj = categorical(ii*ones(size(delta_duration)));
        trialtype = categorical((trialnumber>455)+1) ; % Hill = 1, Valley = 2
        % Create Table for Hill Trials
        temptable = table(subj,trialtype,trialnumber,choicetypecount,...
                          delta_duration,delta_onset,delta_offset,...
                          choseLeft);
        % Remove bad trials (NaNs in the choseLeft field)
        temptable(any(ismissing(temptable),2), :) = [];
        % Update to logical now that NaNs are removed
        temptable.choseLeft = logical(temptable.choseLeft);
        % Add to full table
        TABLE = [TABLE;temptable];

end

%% Initialize Figures
figure(301)
title('Hiils')
ylabel('p(Choosing Left)')
xlim([-4 4]);ylim([0 1])
xlabel('shorter hill <--  Delta Duration of Hill --> longer hill')

figure(302)
title('Valleys')
ylabel('p(Choosing Left)')
xlim([-4 4]);ylim([0 1])
xlabel('shorter valley <--  Delta Duration of Valley --> longer valley')


%% Plot Average across subjects
mean_pChoice_delta_duration_all_subj = nan(2,length(unique_delta_duration));
se_pChoice_delta_duration_all_subj = nan(2,length(unique_delta_duration));
mean_pChoice_delta_onset_all_subj = nan(2,length(unique_delta_duration));
se_pChoice_delta_onset_all_subj = nan(2,length(unique_delta_duration));
mean_pChoice_delta_onset_byDur_all_subj = [];
se_pChoice_delta_onset_byDur_all_subj = [];

for jj = 1:2 % Hill or Valley
    for kk = 1:length(unique_delta_duration)
        mean_pChoice_delta_duration_all_subj(jj,kk) = nanmean(TABLE((TABLE.delta_duration == unique_delta_duration(kk)) & ...
        (TABLE.trialtype == categorical(jj)) & ...
        (TABLE.choicetypecount >= nChoice),:).choseLeft);
        se_pChoice_delta_duration_all_subj(jj,kk) = nanstd(TABLE((TABLE.delta_duration == unique_delta_duration(kk)) & ...
        (TABLE.trialtype == categorical(jj)) & ...
        (TABLE.choicetypecount >= nChoice),:).choseLeft)/sqrt(length(subjs));
    end
    
    figure(300+jj)
    hold on
    plot(unique_delta_duration,mean_pChoice_delta_duration_all_subj(jj,:),'Color',colorvec(jj,:),'LineWidth',3)
    errorbar(unique_delta_duration,mean_pChoice_delta_duration_all_subj(jj,:),se_pChoice_delta_duration_all_subj(jj,:),'k','Linewidth',2,'LineStyle','None')%SE
    for kk = 1:length(unique_delta_onset)
        mean_pChoice_delta_onset_all_subj(jj,kk) = nanmean(TABLE((TABLE.delta_onset == unique_delta_onset(kk)) & ...
        (TABLE.trialtype == categorical(jj)) & ...
        (TABLE.choicetypecount >= nChoice),:).choseLeft);
        se_pChoice_delta_onset_all_subj(jj,kk) = nanstd(TABLE((TABLE.delta_onset == unique_delta_onset(kk)) & ...
        (TABLE.trialtype == categorical(jj)) & ...
        (TABLE.choicetypecount >= nChoice),:).choseLeft)/sqrt(length(subjs));
    end
   
    for nn = 1:length(unique_delta_duration)
        for kk = 1:length(unique_delta_onset)
            mean_pChoice_delta_onset_byDur_all_subj(jj,nn,kk) = nanmean(TABLE((TABLE.delta_onset == unique_delta_onset(kk)) & ...
            (TABLE.trialtype == categorical(jj)) & ...
            (TABLE.delta_duration == unique_delta_duration(nn)) & ...
            (TABLE.choicetypecount >= nChoice),:).choseLeft);

            n_trials_delta_onset_byDur_all_subj(jj,nn,kk) = size(TABLE((TABLE.delta_onset == unique_delta_onset(kk)) & ...
            (TABLE.trialtype == categorical(jj)) & ...
            (TABLE.delta_duration == unique_delta_duration(nn)) & ...
            (TABLE.choicetypecount >= nChoice),:).choseLeft,1);

            se_pChoice_delta_onset_byDur_all_subj(jj,nn,kk) = nanstd(TABLE((TABLE.delta_onset == unique_delta_onset(kk)) & ...
            (TABLE.trialtype == categorical(jj)) & ...
            (TABLE.delta_duration == unique_delta_duration(nn)) & ...
            (TABLE.choicetypecount >= nChoice),:).choseLeft)/sqrt(length(subjs));
        end
       end
end

plot_colors = cool(length(unique_delta_duration));
for jj = 1:2 % Hill or Valley
%     subplot(2,1,jj)
    figure(302+jj)
    hold on
    for kk = 1:length(unique_delta_duration)
%         ind = ~isnan(squeeze(pChoice_delta_onset_byDur_all_subj(jj,kk,:)));
        ind = ~isnan(squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,:))); %NEW EXP?
        plot(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'color',plot_colors(kk,:),'LineWidth',2)
        errorbar(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),squeeze(se_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'k','Linewidth',2,'LineStyle','None')%SE
        %         plot(unique_delta_onset(ind),squeeze(pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'.','color',plot_colors(kk,:),'MarkerSize',10)
    end
    title(session_strings{jj})
    xlabel('Difference in Onset Time (sec)')
    ylabel('P(Choosing Left)')
    legend('Ref 4s shorter','','Ref 2s shorter','','Ref equal dur','','Ref 2s longer','','Ref 4s longer','')
end


%% Plot Model Results

session_strings = {'Hill','Valley'};
unique_delta_duration = -4:2:4; 
unique_delta_onset = -10:2:10;
pChoice_delta_duration = nan(length(subjs),2,length(unique_delta_duration));
pChoice_delta_onset = nan(length(subjs),2,length(unique_delta_onset));


%% Run Linear Model

choice_table = [];
for jj = 1:length(subjs)
    for kk = 1:length(Mchoice{jj}.choice)
        if ~isnan(Mchoice{jj}.choice(kk))
            choice_table = [choice_table;
                            %subj   hill/vall                       diff duration                                               diff onset                                             chose left ...
                            jj,     Mchoice{jj}.trialnumber(kk) <= 455, Mchoice{jj}.Ltopduration(kk) - Mchoice{jj}.Rtopduration(kk),Mchoice{jj}.Ltoponset(kk) - Mchoice{jj}.Rtoponset(kk), -1*Mchoice{jj}.choice(kk)+2];
        end
    end
end

% Define Variable Types and Name Variables
choice_table = array2table(choice_table,'VariableNames',{'subj','hill_or_valley','delta_duration','delta_onset','chose_left'});
choice_table.subj = categorical(choice_table.subj);
choice_table.hill_or_valley = categorical(choice_table.hill_or_valley);
choice_table.chose_left = logical(choice_table.chose_left);

% Fit Model
mdl_choice = fitglme(choice_table,'chose_left ~ hill_or_valley*delta_duration*delta_onset + (1|subj)','Distribution','Binomial','Link','Logit')

%% Plot it all together

% Start with Difference in duration
%Hill
model_table = table();
model_table.delta_duration = [-4:0.5:4]';
model_table.subj = categorical(zeros(size(model_table.delta_duration)));
model_table.hill_or_valley = categorical(ones(size(model_table.delta_duration))); % Hills First
model_table.delta_onset = zeros(size(model_table.delta_duration));

figure(301)
[y_pred_delta_dur_hill] = predict(mdl_choice,model_table);
plot(model_table.delta_duration,y_pred_delta_dur_hill,':','Color','r','LineWidth',3)

% Now Valleys
figure(302)
model_table.hill_or_valley = categorical(zeros(size(model_table.delta_duration))); % Valleys
[y_pred_delta_dur_valley] = predict(mdl_choice,model_table);
plot(model_table.delta_duration,y_pred_delta_dur_valley,':','Color','r','LineWidth',3)

% Now Difference in Onset Time
% possible_combos_hill = [4*ones(1,9),2*ones(1,10),zeros(1,10),-2*ones(1,10),-4*ones(1,9);
%                         -10:2:6,-10:2:8,-10:2:-2,2:2:10,-8:2:10,-6:2:10]';

possible_combos = [4*ones(1,9),2*ones(1,10),zeros(1,10),-2*ones(1,10),-4*ones(1,9);
                        -10:2:6,-10:2:8,-10:2:-2,2:2:10,-8:2:10,-6:2:10]';

model_table = table();
model_table.delta_duration = possible_combos(:,1);
model_table.delta_onset = possible_combos(:,2);
model_table.subj = categorical(zeros(size(possible_combos,1),1));
model_table.hill_or_valley = categorical(ones(size(possible_combos,1),1)); % Hills First
[y_pred_delta_onset_hill] = predict(mdl_choice,model_table);

plot_colors = cool(length(unique(model_table.delta_duration)));

% Plot By Diff in Duration
delta_durations = -4:2:4;
% gcf;subplot(2,1,1);hold on;
figure(303)
for ii = 1:length(delta_durations)
    plot(possible_combos(possible_combos(:,1) == delta_durations(ii),2),y_pred_delta_onset_hill(possible_combos(:,1) == delta_durations(ii)),':','Color',plot_colors(ii,:),'LineWidth',3)
end
title('Hills')

% Now Valleys
model_table.hill_or_valley = categorical(zeros(size(possible_combos,1),1)); % Hills First
[y_pred_delta_onset_hill] = predict(mdl_choice,model_table);

% Plot By Diff in Duration
delta_durations = -4:2:4;
% subplot(2,1,2);hold on;
figure(304)
for ii = 1:length(delta_durations)
    plot(possible_combos(possible_combos(:,1) == delta_durations(ii),2),y_pred_delta_onset_hill(possible_combos(:,1) == delta_durations(ii)),':','Color',plot_colors(ii,:),'LineWidth',3)
end
title('Valleys')



%% Plot All Subjects, for Hill / Valley When Durations are equal

possible_combos = [4*ones(1,9),2*ones(1,10),zeros(1,10),-2*ones(1,10),-4*ones(1,9);
                        -10:2:6,-10:2:8,-10:2:-2,2:2:10,-8:2:10,-6:2:10]';

model_table = table();
model_table.delta_duration = possible_combos(:,1);
model_table.delta_onset = possible_combos(:,2);
model_table.subj = categorical(zeros(size(possible_combos,1),1));
model_table.hill_or_valley = categorical(ones(size(possible_combos,1),1)); % Hills First
[y_pred_delta_onset_hill] = predict(mdl_choice,model_table);

plot_colors = cool(length(unique(model_table.delta_duration)));

pChoices = nan(length(subjs),2,length(unique_delta_duration),length(unique_delta_onset));
for jj = 1:2 % Hill or Valley
    for kk = 1:length(unique_delta_duration)
        for mm = 1:length(unique_delta_onset)
%         delta_onsets_temp = possible_combos(possible_combos(:,1) == delta_durations(3),2);
            for nn = 1:length(subjs)
%             for mm = 1:length(delta_onsets_temp)
                pChoices(nn,jj,kk,mm) = nanmean(TABLE((TABLE.delta_duration == unique_delta_duration(kk)) & ...
                (TABLE.trialtype == categorical(jj)) & ...
                (TABLE.subj == categorical(nn)) & ...
                (TABLE.delta_onset == unique_delta_onset(mm)) & ...
                (TABLE.choicetypecount >= nChoice),:).choseLeft);
            end
        end
    end
end


% Hills
figure(305)
hold on
title('Hills')
% All Individuals
for nn = 1:length(subjs)
    ind = ~isnan(squeeze(pChoices(nn,1,3,:)));
    plot(unique_delta_onset(ind),squeeze(pChoices(nn,1,3,ind)),'LineWidth',2,'Color',[0,0,0,0.05])    
    %     plot(possible_combos(possible_combos(:,1) == delta_durations(3),2),rmmissing(squeeze(pChoices(nn,1,3,:))),'LineWidth',2,'Color',[0,0,0,0.05])
end
% Mean Subj Data
ind = ~isnan(squeeze(mean_pChoice_delta_onset_byDur_all_subj(1,3,:))); %NEW EXP?
plot(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(1,3,ind)),'color',plot_colors(3,:),'LineWidth',2)
errorbar(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(1,3,ind)),squeeze(se_pChoice_delta_onset_byDur_all_subj(1,3,ind)),'k','Linewidth',2,'LineStyle','None')%SE
% Model Prediction
plot(possible_combos(possible_combos(:,1) == delta_durations(3),2),y_pred_delta_onset_hill(possible_combos(:,1) == delta_durations(3)),':','Color',plot_colors(3,:),'LineWidth',3)


% Valleys
model_table.hill_or_valley = categorical(zeros(size(possible_combos,1),1)); 
[y_pred_delta_onset_valley] = predict(mdl_choice,model_table);

figure(306)
hold on
title('Valleys')
% All Individuals
for nn = 1:length(subjs)
    ind = ~isnan(squeeze(pChoices(nn,2,3,:)));
    plot(unique_delta_onset(ind),squeeze(pChoices(nn,2,3,ind)),'LineWidth',2,'Color',[0,0,0,0.05])    
%     plot(possible_combos(possible_combos(:,1) == delta_durations(3),2),rmmissing(squeeze(pChoices(nn,2,3,:))),'LineWidth',2,'Color',[0,0,0,0.05])
end
% Mean Subj Data
ind = ~isnan(squeeze(mean_pChoice_delta_onset_byDur_all_subj(2,3,:))); %NEW EXP?
plot(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(2,3,ind)),'color',plot_colors(3,:),'LineWidth',2)
errorbar(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(2,3,ind)),squeeze(se_pChoice_delta_onset_byDur_all_subj(2,3,ind)),'k','Linewidth',2,'LineStyle','None')%SE
% Model Prediction
plot(possible_combos(possible_combos(:,1) == delta_durations(3),2),y_pred_delta_onset_valley(possible_combos(:,1) == delta_durations(3)),':','Color',plot_colors(3,:),'LineWidth',3)



for jj = 1:2 % Hill or Valley
%     subplot(2,1,jj)
    figure(302+jj)
    hold on
    for kk = 1:length(unique_delta_duration)
%         ind = ~isnan(squeeze(pChoice_delta_onset_byDur_all_subj(jj,kk,:)));
        ind = ~isnan(squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,:))); %NEW EXP?
        plot(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'color',plot_colors(kk,:),'LineWidth',2)
        errorbar(unique_delta_onset(ind),squeeze(mean_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),squeeze(se_pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'k','Linewidth',2,'LineStyle','None')%SE
        %         plot(unique_delta_onset(ind),squeeze(pChoice_delta_onset_byDur_all_subj(jj,kk,ind)),'.','color',plot_colors(kk,:),'MarkerSize',10)
    end
    title(session_strings{jj})
    xlabel('Difference in Onset Time (sec)')
    ylabel('P(Choosing Left)')
    legend('Ref 4s shorter','','Ref 2s shorter','','Ref equal dur','','Ref 2s longer','','Ref 4s longer','')
end
