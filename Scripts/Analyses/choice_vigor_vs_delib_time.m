% test correlation of delib time, choice vigor (and pre- or pro-
% crastination)

%% Process Subject Fits
time_physical = squeeze(sol(1,4,:));
time_tasks = squeeze(sol(2,4,:));
time_physical_cat_ispro = (time_physical>0);
time_tasks_cat_ispro = (time_tasks>0);

% include_outlier_trials = 0; 
% include_outlier_subjs = 1;
% norm_to_mvc = 1;
norm_to_mvc = 0;


selectedmodel = 4*ones(length(subjs),1);

% Based on Best Fit Model Structure for Majority of Subjects
% selectedmodel = modelstocompare(find(bmsout.pxp==max(bmsout.pxp)))*ones(length(subjs),1);
% OR Based on Best Fit Model for that Subject
% selectedmodel = bestmodel;

sol_reduced = squeeze(sol); % should be (model, subj, param)

choiceforcerate_table = [];
% For Each Subject
for jj = 1:length(subjs)

%     gamma_physical = sol_reduced(selectedmodel(jj),jj,1);
    gamma_physical = sol_reduced(1,selectedmodel(jj),jj);
%     gamma_task = sol_reduced(selectedmodel(jj),jj,2);
    gamma_task = sol_reduced(2,selectedmodel(jj),jj);
%     beta = sol_reduced(selectedmodel(jj),jj,3);
    beta = sol_reduced(3,selectedmodel(jj),jj);
%     rel_weight = sol_reduced(selectedmodel(jj),jj,4); 
    if selectedmodel(jj) <=2
        rel_weight = 1;
    else
        rel_weight = sol_reduced(4,selectedmodel(jj),jj); 
    end

    Mchoice{jj}.subj_u_L = nan(Mchoice{jj}.trials,1);
    Mchoice{jj}.subj_u_R = nan(Mchoice{jj}.trials,1);
    Mchoice{jj}.delta_utility = nan(Mchoice{jj}.trials,1);
    Mchoice{jj}.p_predicted = nan(Mchoice{jj}.trials,1);

    % Calculate 
    for ii = 1:Mchoice{jj}.trials
        Mchoice{jj}.subj_u_L(ii) = rel_weight*subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice{jj}.Ltrialduration(ii),Mchoice{jj}.Lrampuptime2base(ii),Mchoice{jj}.Lrampuptime2top(ii),...
                                Mchoice{jj}.Lrampdowntime2base(ii),Mchoice{jj}.Lrampdowntime2end(ii),Mchoice{jj}.Lbaselevel(ii),Mchoice{jj}.Ltoplevel(ii),...
                                Mchoice{jj}.Ltopduration(ii),Mchoice{jj}.Ltoponset(ii),discounttype) + ...
           (1-rel_weight)*subjectiveutility_task(gamma_task,reward,Mchoice{jj}.Ltrialduration(ii),Mchoice{jj}.Lrampuptime2base(ii),Mchoice{jj}.Lrampuptime2top(ii),...
                                Mchoice{jj}.Lrampdowntime2base(ii),Mchoice{jj}.Lrampdowntime2end(ii),Mchoice{jj}.Lbaselevel(ii),Mchoice{jj}.Ltoplevel(ii),...
                                Mchoice{jj}.Ltopduration(ii),Mchoice{jj}.Ltoponset(ii),taskmodeltype,discounttype);    
        Mchoice{jj}.subj_u_R(ii) = rel_weight*subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice{jj}.Rtrialduration(ii),Mchoice{jj}.Rrampuptime2base(ii),Mchoice{jj}.Rrampuptime2top(ii),...
                                Mchoice{jj}.Rrampdowntime2base(ii),Mchoice{jj}.Rrampdowntime2end(ii),Mchoice{jj}.Rbaselevel(ii),Mchoice{jj}.Rtoplevel(ii),...
                                Mchoice{jj}.Rtopduration(ii),Mchoice{jj}.Rtoponset(ii),discounttype) + ...
           (1-rel_weight)*subjectiveutility_task(gamma_task,reward,Mchoice{jj}.Rtrialduration(ii),Mchoice{jj}.Rrampuptime2base(ii),Mchoice{jj}.Rrampuptime2top(ii),...
                                Mchoice{jj}.Rrampdowntime2base(ii),Mchoice{jj}.Rrampdowntime2end(ii),Mchoice{jj}.Rbaselevel(ii),Mchoice{jj}.Rtoplevel(ii),...
                                Mchoice{jj}.Rtopduration(ii),Mchoice{jj}.Rtoponset(ii),taskmodeltype,discounttype); 
        if ~isnan(Mchoice{jj}.choice(ii))
            % Probability of Choosing Left Choice
            Mchoice{jj}.p_predicted(ii) = ( 1 + exp(-1*beta*(Mchoice{jj}.subj_u_L(ii) - Mchoice{jj}.subj_u_R(ii))) ).^-1;
            if Mchoice{jj}.p_predicted(ii) == 0
                Mchoice{jj}.p_predicted(ii) = eps;
            elseif Mchoice{jj}.p_predicted(ii) == 1
                Mchoice{jj}.p_predicted(ii) = 1-eps;
            end
            % Delta Utility
            Mchoice{jj}.delta_utility(ii) = Mchoice{jj}.subj_u_L(ii) - Mchoice{jj}.subj_u_R(ii);
        else
            Mchoice{jj}.p_predicted(ii) = NaN;
            Mchoice{jj}.delta_utility(ii) = NaN;
        end
    end
    Mchoice{jj}.delta_utility_normed = Mchoice{jj}.delta_utility/max(Mchoice{jj}.delta_utility);
    Mchoice{jj}.abs_delta_utility_minimaxed = (abs(Mchoice{jj}.delta_utility) - min(abs(Mchoice{jj}.delta_utility)))/(max(abs(Mchoice{jj}.delta_utility))-min(abs(Mchoice{jj}.delta_utility)));
end

%% Prepare Data Table
choice_table = [];
subject_outliers = [];
% Create new array without outlier trials within subject, Delib
if include_outlier_trials == 0
    for jj = 1:length(subjs)
        Mchoice{jj}.delibtime_nooutliers = Mchoice{jj}.delibtime;
%             Mchoice{jj}.delibtime_nooutliers(isoutlier(Mchoice{jj}.delibtime,'percentiles',[1 99])) = NaN;
        Mchoice{jj}.delibtime_nooutliers(isoutlier(Mchoice{jj}.delibtime,'percentiles',[2.5 97.5])) = NaN;
        Mchoice{jj}.delibtime_nooutliers(Mchoice{jj}.delibtime_nooutliers>5) = NaN;

        mean_delib(jj) = nanmean(Mchoice{jj}.delibtime_nooutliers);
        std_delib(jj) = nanstd(Mchoice{jj}.delibtime_nooutliers);
    end
else 
    for jj = 1:length(subjs)
        mean_delib(jj) = nanmean(Mchoice{jj}.delibtime);
        std_delib(jj) = nanstd(Mchoice{jj}.delibtime);
    end
end
% Vigor
% Create new array without outlier trials within subject
if include_outlier_trials == 0
    for jj = 1:length(subjs)
        Mchoice{jj}.max_dxpx_dt_nooutliers = Mchoice{jj}.max_dxpx_dt;
%             Mchoice{jj}.max_dxpx_dt_nooutliers(isoutlier(Mchoice{jj}.max_dxpx_dt,'percentiles',[1 99])) = NaN;
        Mchoice{jj}.max_dxpx_dt_nooutliers(isoutlier(Mchoice{jj}.max_dxpx_dt,'percentiles',[2.5 97.5])) = NaN;
        mean_vigor(jj) = nanmean(Mchoice{jj}.max_dxpx_dt_nooutliers);
        std_vigor(jj) = nanstd(Mchoice{jj}.max_dxpx_dt_nooutliers);

    end
else 
    for jj = 1:length(subjs)
        mean_vigor(jj) = nanmean(Mchoice{jj}.max_dxpx_dt);
        std_vigor(jj) = nanstd(Mchoice{jj}.max_dxpx_dt);
    end
end

% Determine if a subject is an outlier based on deliberation time
subject_outliers = (isoutlier(mean_delib) | isoutlier(mean_vigor));
if include_outlier_subjs == 1
    nsubjs = length(subjs);
else
    nsubjs = sum(subject_outliers==0);
end


for jj = 1:length(subjs)
    % Z-Score Data
%         Mchoice{jj}.delib_time_zscored = (Mchoice{jj}.delibtime - mean_delib(jj))./std_delib(jj);
    for kk = 1:Mchoice{jj}.trials
        % Move this to original processing TODO
        if kk <= 455
            session = 1; %A
        else
            session = 2; %B
        end
        if (subject_outliers(jj) == 0 && include_outlier_subjs == 0) || include_outlier_subjs == 1
            if include_outlier_trials == 0 && ~isnan(Mchoice{jj}.delibtime_nooutliers(kk)) && ~isnan(Mchoice{jj}.delibtime_nooutliers(kk))
                choice_table = [choice_table;
                                jj,time_physical(jj),time_tasks(jj),time_physical_cat_ispro(jj),time_tasks_cat_ispro(jj),...
                                abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime_nooutliers(kk),Mchoice{jj}.max_dxpx_dt_nooutliers(kk)];
            else
               choice_table = [choice_table;
                                jj,time_physical(jj),time_tasks(jj),time_physical_cat_ispro(jj),time_tasks_cat_ispro(jj),...
                                abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime(kk),Mchoice{jj}.max_dxpx_dt(kk)];
            end
        end
    end
end
choice_table = array2table(choice_table,'VariableNames',{'subj','gamma_physical','gamma_task','physical_ispro','task_ispro',...
                                                                    'abs_dUtil_normed','delib_time','choice_vigor'});
choice_table.subj = categorical(choice_table.subj);
choice_table.physical_ispro = logical(choice_table.physical_ispro);
choice_table.task_ispro = logical(choice_table.task_ispro);

if norm_to_mvc
    choice_table.choice_vigor = (choice_table.choice_vigor)./824;
end

% delib_and_vigor = fitlme(choice_table,'choice_vigor ~ gamma_physical * gamma_task + (1|subj)')
% delib_and_vigor = fitlme(choice_table,'delib_time ~ gamma_physical * gamma_task + (1|subj)')


% delib_and_vigor = fitlme(choice_table,'choice_vigor ~ delib_time + (1|subj)')
delib_and_vigor = fitlme(choice_table,'delib_time ~ choice_vigor + (1|subj)')


% figure
% hold on
% plot(choice_table.delib_time,choice_table.choice_vigor,'x')
% plt_x = min(choice_table.delib_time):2:max(choice_table.delib_time);
% plot(plt_x,delib_and_vigor.Coefficients(1,2).Estimate + delib_and_vigor.Coefficients(2,2).Estimate*plt_x,'k','LineWidth',2)
%

bins = 0:2:12; % if %MVC
bins = bins*824;
[Y_delib_byvigor,~] = discretize(choice_table.choice_vigor,bins);
[meanVal_vigor, maxVal_vigor, stdVal_vigor, medVal_vigor] = grpstats(choice_table.delib_time,Y_delib_byvigor,{@mean, @max, @std, @median});
[x_meanVal_vigor,x_maxVal_vigor, x_stdVal_vigor, x_medVal_vigor] = grpstats(choice_table.choice_vigor,Y_delib_byvigor,{@mean, @max, @std, @median});
figure(503)
hold on
errorbar(x_meanVal_vigor,meanVal_vigor,stdVal_vigor/sqrt(nsubjs),stdVal_vigor/sqrt(nsubjs),...
        x_stdVal_vigor./sqrt(nsubjs),x_stdVal_vigor./sqrt(nsubjs),'k','Linewidth',2,'LineStyle','None')%SE

tblnew = table();
tblnew.choice_vigor = bins';
tblnew.subj = 0*ones(size(tblnew.choice_vigor));
tblnew.subj = categorical(tblnew.subj);
[test,test_CI] = predict(delib_and_vigor,tblnew);
x_CI = [bins,bins(end:-1:1)];
y_CI = [test_CI(:,1)',test_CI(end:-1:1,2)'];
plot(bins,test,'Color',[150,100,150]./255,'LineWidth',3);
p = fill(x_CI,y_CI,'red');p.FaceColor = [150,100,150]./255;p.EdgeColor = 'None';p.FaceAlpha = 0.5;
ylim([0 12*824])
xlim([0 12*824])
ylim([0 1.2])
xlabel('Maximum Cursor Velocity (px/sec)')
ylabel('Deliberation Time (s)')
title('(C) - Deliberation Time vs. Choice Vigor')


