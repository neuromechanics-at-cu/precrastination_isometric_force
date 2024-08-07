
%% Investigate Realization Trial Performance

table = [];

for jj = 1:length(subjs)
    for kk = 1:Mchoice{jj}.trials
        if Mchoice{jj}.realization(kk)

            if Mchoice{jj}.choice(kk) == 1
                duration_temp = Mchoice{jj}.Ltopduration(kk);
                onset_temp = Mchoice{jj}.Ltoponset(kk);
            elseif Mchoice{jj}.choice(kk) == 2
                duration_temp = Mchoice{jj}.Rtopduration(kk);
                onset_temp = Mchoice{jj}.Rtoponset(kk);
            end

            table = [table;
                    jj,kk<=455,Mchoice{jj}.trialnumber_bysession(kk), Mchoice{jj}.trialnumber_bysession(kk)<=227,duration_temp,onset_temp,Mchoice{jj}.avg_force_base1(kk),Mchoice{jj}.std_force_base1(kk),Mchoice{jj}.avg_force_hillvalley(kk),Mchoice{jj}.std_force_hillvalley(kk),...
                    Mchoice{jj}.avg_force_base2(kk),Mchoice{jj}.std_force_base2(kk),Mchoice{jj}.overshoot_ramp1(kk),Mchoice{jj}.peakrate_mvc_ramp1(kk),Mchoice{jj}.overshoot_ramp2(kk),Mchoice{jj}.peakrate_mvc_ramp2(kk)];
        end
    end
end
            
%% Adjust Table
% table = array2table(table,'VariableNames',{'subj','hill_or_valley','trial_num','duration','onset','avg_base_1','std_base_1','avg_hillvall','std_hillvall','avg_base_2','std_base_2','overshoot_ramp_1','peak_rate_ramp_1','overshoot_ramp_2','peak_rate_ramp_2'});
table = array2table(table,'VariableNames',{'subj','hill_or_valley','trial_num','first_half','duration','onset','avg_base_1','std_base_1','avg_hillvall','std_hillvall','avg_base_2','std_base_2','overshoot_ramp_1','peak_rate_ramp_1','overshoot_ramp_2','peak_rate_ramp_2'});
table.subj = categorical(table.subj);
table.hill_or_valley = categorical(table.hill_or_valley);

fitlme(table,'avg_base_1 ~ hill_or_valley*trial_num*duration*onset + (1|subj)')
fitlme(table,'avg_base_2 ~ hill_or_valley*trial_num*duration*onset + (1|subj)')
fitlme(table,'avg_hillvall ~ hill_or_valley*trial_num*duration*onset + (1|subj)')