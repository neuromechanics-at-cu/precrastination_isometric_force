% taskmodeltype = 0; 
% discounttype = 2;

%% Calculate Subjective Utility Based on Model Fits
AICs = [];
BICs = [];
Ftests = [];
AICs2 = [];

% include_outlier_trials = 0; 
% include_outlier_subjs = 1;
% 
% 
% modelstocompare = 1:4;

for nn = 1:length(modelstocompare)
    selectedmodel = modelstocompare(nn)*ones(length(subjs),1);

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
    choice_time_table = [];
    subject_outliers = [];
    % Create new array without outlier trials within subject
    if include_outlier_trials == 0
        for jj = 1:length(subjs)
            Mchoice{jj}.delibtime_nooutliers = Mchoice{jj}.delibtime;
            Mchoice{jj}.delibtime_nooutliers(isoutlier(Mchoice{jj}.delibtime,'percentiles',[2.5 97.5])) = NaN;
            Mchoice{jj}.delibtime_nooutliers((Mchoice{jj}.delibtime_nooutliers>4)) = NaN;
            mean_delib(jj) = nanmean(Mchoice{jj}.delibtime_nooutliers);
            std_delib(jj) = nanstd(Mchoice{jj}.delibtime_nooutliers);
        end
    else 
        for jj = 1:length(subjs)
            mean_delib(jj) = nanmean(Mchoice{jj}.delibtime);
            std_delib(jj) = nanstd(Mchoice{jj}.delibtime);
        end
    end
    % Determine if a subject is an outlier
    subject_outliers = isoutlier(mean_delib);
    if include_outlier_subjs == 1
        nsubjs = length(subjs);
    else
        nsubjs = sum(subject_outliers==0);
    end
    
    for jj = 1:length(subjs)
        % Z-Score Data
        Mchoice{jj}.delib_time_zscored = (Mchoice{jj}.delibtime - mean_delib(jj))./std_delib(jj);
        Mchoice{jj}.delib_time_zscored_nooutliers = (Mchoice{jj}.delibtime_nooutliers - mean_delib(jj))./std_delib(jj);

        if (subject_outliers(jj) == 0 && include_outlier_subjs == 0) || include_outlier_subjs == 1
            if include_outlier_trials == 0
                for kk = 1:length(Mchoice{jj}.delibtime_nooutliers)
                    if ~isnan(Mchoice{jj}.delibtime_nooutliers(kk))
                        choice_time_table = [choice_time_table;
%                                          jj,abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime_nooutliers(kk),Mchoice{jj}.delibtime_nooutliers_normed(kk),Mchoice{jj}.delib_time_zscored_nooutliers(kk)];
                                         jj,Mchoice{jj}.trialnumber_bysession(kk),abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime_nooutliers(kk),Mchoice{jj}.delibtime_nooutliers_normed(kk),Mchoice{jj}.delib_time_zscored_nooutliers(kk)];
                
                    end
                end
            else 
                for kk = 1:length(Mchoice{jj}.delibtime)
                    choice_time_table = [choice_time_table;
%                                      jj,abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime(kk),Mchoice{jj}.delibtime_normed(kk),Mchoice{jj}.delib_time_zscored(kk)];
                                     jj,Mchoice{jj}.trialnumber_bysession(kk),abs(Mchoice{jj}.abs_delta_utility_minimaxed(kk)),Mchoice{jj}.delibtime(kk),Mchoice{jj}.delibtime_normed(kk),Mchoice{jj}.delib_time_zscored(kk)];
    
                end    
            end
        end
    end
%     choice_time_table = array2table(choice_time_table,'VariableNames',{'subj','abs_dUtil_normed','delib_time','normed_delib_time','zscored_delib_time'});
    choice_time_table = array2table(choice_time_table,'VariableNames',{'subj','trial_num','abs_dUtil_normed','delib_time','normed_delib_time','zscored_delib_time'});
    choice_time_table.subj = categorical(choice_time_table.subj);

    glme_mdl = fitlme(choice_time_table,'delib_time ~ abs_dUtil_normed + (1|subj)')

    bin_midpts = 0:.2:1;
    % Bin data and plot for visualization
    [Y,E] = discretize(choice_time_table.abs_dUtil_normed,bin_midpts);

    
    % Plot binned results

    [meanVal, maxVal, stdVal] = grpstats(choice_time_table.delib_time,Y,{@mean, @max, @std}); 
    [x_meanVal,x_maxVal, x_stdVal] = grpstats(choice_time_table.abs_dUtil_normed,Y,{@mean, @max, @std}); 
    

%     bin_midpts = 0:.1:1;
% FEWER BINS
%     bin_midpts = 0.1:.2:0.9;
    domain = 0:0.05:1;


    
    % Plot Along the Way
    figure(510+nn)
    hold on

    
    tblnew = table();
    tblnew.abs_dUtil_normed = domain';
    tblnew.subj = 0*ones(size(tblnew.abs_dUtil_normed));
    tblnew.subj = categorical(tblnew.subj);
    [test,test_CI] = predict(glme_mdl,tblnew);
    plot(domain,test,'Color',[150,100,150]./255,'LineWidth',3);
    x_CI = [domain,domain(end:-1:1)];
    y_CI = [test_CI(:,1)',test_CI(end:-1:1,2)'];
    p = fill(x_CI,y_CI,'red');p.FaceColor = [150,100,150]./255;p.EdgeColor = 'None';p.FaceAlpha = 0.5;

    errorbar(x_meanVal,meanVal,stdVal/sqrt(nsubjs),stdVal/sqrt(nsubjs),...
        x_stdVal./sqrt(nsubjs),x_stdVal./sqrt(nsubjs),'k','Linewidth',2,'LineStyle','None')%SE

    xlabel('Normalized Difference in Subjective Utility')
    ylabel('Choice Deliberation Time (s)')
    title('(A) Deliberation Time vs. Difference in Utility')



    AICs = [AICs,glme_mdl.ModelCriterion.AIC];
    BICs = [BICs,glme_mdl.ModelCriterion.BIC];
    Ftests = [Ftests,coefTest(glme_mdl)];
    MDLS{nn} = glme_mdl;
end

%% Print AICS
AICs
%% Print Model Parameters


%% Dont Plot Losing Models
[~,ind_winner] = min(AICs);

if 1 
    for nn = 1:length(modelstocompare)
        if nn ~= ind_winner
            close(figure(510+nn))
        end
    end
end

