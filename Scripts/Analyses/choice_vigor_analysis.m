
%% Load in Model Data


% discounttype = 2; %% DELETE THIS TOO TODO
% taskmodeltype = 0; %% DELETE THIS TOO TODO


addpath('../Utility Model Fitting/')
addpath('../Supporting Functions')


%% Calculate Subjective Utility Based on Model Fits
AICs = [];
BICs = [];
Ftests = [];


% include_outlier_trials = 0; 
% include_outlier_subjs = 1;
norm_to_mvc = 1;

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

        %% Compare Vigor Correlations
        Mchoice{jj}.subj_u_L_normed = -1*Mchoice{jj}.subj_u_L/min([Mchoice{jj}.subj_u_L',Mchoice{jj}.subj_u_R']);
        Mchoice{jj}.subj_u_R_normed = -1*Mchoice{jj}.subj_u_R/min([Mchoice{jj}.subj_u_L',Mchoice{jj}.subj_u_R']);

        %% Transform Max Rate to %MVC
%         Mchoice{jj}.max_dxpx_dt = Mchoice{jj}.max_dxpx_dt/824;
    end

    subject_outliers = [];
    % Create new array without outlier trials within subject
    if include_outlier_trials == 0
        for jj = 1:length(subjs)
            Mchoice{jj}.max_dxpx_dt_nooutliers = Mchoice{jj}.max_dxpx_dt;
            Mchoice{jj}.max_dxpx_dt_nooutliers(isoutlier(Mchoice{jj}.max_dxpx_dt,'percentiles',[2.5 97.5])) = NaN;
            mean_vigor(jj) = nanmean(Mchoice{jj}.max_dxpx_dt_nooutliers);
            std_vigor(jj) = nanstd(Mchoice{jj}.max_dxpx_dt_nooutliers);
            Mchoice{jj}.max_dxpx_dt_zscored_nooutliers = (Mchoice{jj}.max_dxpx_dt_nooutliers - mean_vigor(jj))/std_vigor(jj);
            Mchoice{jj}.max_dxpx_dt_normed_nooutliers = Mchoice{jj}.max_dxpx_dt_nooutliers - mean_vigor(jj);
        end
    else 
        for jj = 1:length(subjs)
            mean_vigor(jj) = nanmean(Mchoice{jj}.max_dxpx_dt);
            std_vigor(jj) = nanstd(Mchoice{jj}.max_dxpx_dt);
            Mchoice{jj}.max_dxpx_dt_zscored = (Mchoice{jj}.max_dxpx_dt - mean_vigor(jj))/std_vigor(jj);
            Mchoice{jj}.max_dxpx_dt_normed = Mchoice{jj}.max_dxpx_dt - mean_vigor(jj);
        end
    end
    % Determine if a subject is an outlier
    subject_outliers = isoutlier(mean_vigor);
    if include_outlier_subjs == 1
        nsubjs = length(subjs);
    else
        nsubjs = sum(subject_outliers==0);
    end
    for jj = 1:length(subjs)

        if (subject_outliers(jj) == 0 && include_outlier_subjs == 0) || include_outlier_subjs == 1
    %         Mchoice{jj}.delibtime_normedpersubj = Mchoice{jj}.delibtime_nooutliers - mean(Mchoice{jj}.delibtime_nooutliers);
            if include_outlier_trials == 0
                for kk = 1:length(Mchoice{jj}.max_dxpx_dt_nooutliers)
                    if ~isnan(Mchoice{jj}.delibtime_nooutliers(kk))
%                         choiceforcerate_table = [choiceforcerate_table;
%                                              jj,Mchoice{jj}.abs_delta_utility_minimaxed(kk),...
%                                              Mchoice{jj}.max_dxpx_dt_nooutliers(kk),log(Mchoice{jj}.max_dxpx_dt_nooutliers(kk)),...
%                                              Mchoice{jj}.max_dxpx_dt_normed_nooutliers(kk),Mchoice{jj}.max_dxpx_dt_zscored_nooutliers(kk)];
                        choiceforcerate_table = [choiceforcerate_table;
                                             jj,Mchoice{jj}.trialnumber_bysession(kk),Mchoice{jj}.abs_delta_utility_minimaxed(kk),...
                                             Mchoice{jj}.max_dxpx_dt_nooutliers(kk),log(Mchoice{jj}.max_dxpx_dt_nooutliers(kk)),...
                                             Mchoice{jj}.max_dxpx_dt_normed_nooutliers(kk),Mchoice{jj}.max_dxpx_dt_zscored_nooutliers(kk)];                    
                    end
                end
            else 
                for kk = 1:length(Mchoice{jj}.delibtime)
%                     choiceforcerate_table = [choiceforcerate_table;
%                                              jj,Mchoice{jj}.abs_delta_utility_minimaxed(kk),...
%                                              Mchoice{jj}.max_dxpx_dt(kk),log(Mchoice{jj}.max_dxpx_dt(kk)),...
%                                              Mchoice{jj}.max_dxpx_dt_normed(kk),Mchoice{jj}.max_dxpx_dt_zscored(kk)];    
                    choiceforcerate_table = [choiceforcerate_table;
                                             jj,Mchoice{jj}.trialnumber_bysession(kk),Mchoice{jj}.abs_delta_utility_minimaxed(kk),...
                                             Mchoice{jj}.max_dxpx_dt(kk),log(Mchoice{jj}.max_dxpx_dt(kk)),...
                                             Mchoice{jj}.max_dxpx_dt_normed(kk),Mchoice{jj}.max_dxpx_dt_zscored(kk)];                    
                end    
            end
        end

    end
%     choiceforcerate_table = array2table(choiceforcerate_table,'VariableNames',{'subj','abs_dUtil_normed','abs_max_rate','log_abs_max_rate','normed_max_rate','zscored_max_rate'});
    choiceforcerate_table = array2table(choiceforcerate_table,'VariableNames',{'subj','trial_num','abs_dUtil_normed','abs_max_rate','log_abs_max_rate','normed_max_rate','zscored_max_rate'});
  
    choiceforcerate_table.subj = categorical(choiceforcerate_table.subj);

    % Normalize to %MVC/sec
    if norm_to_mvc
        choiceforcerate_table.abs_max_rate = (choiceforcerate_table.abs_max_rate)./824;
    end
    mdl_delta_util = fitlme(choiceforcerate_table,'abs_max_rate ~ abs_dUtil_normed + (1|subj)') % USE THIS ONE
%     bin_midpts = 0:.1:1;
    bin_midpts = 0:.2:1;

    plot_pts = [0,0.5,1];
%     [Y_vigor,~] = discretize(choiceforcerate_table.abs_dUtil_normed,0:.1:1,bin_midpts);
    [Y_vigor,~] = discretize(choiceforcerate_table.abs_dUtil_normed,bin_midpts);
%     [meanVal_vigor, maxVal_vigor, stdVal_vigor, medVal_vigor] = grpstats(choiceforcerate_table.zscored_max_rate,Y_vigor,{@mean, @max, @std, @median}); 
    [meanVal_vigor, maxVal_vigor, stdVal_vigor, medVal_vigor] = grpstats(choiceforcerate_table.abs_max_rate,Y_vigor,{@mean, @max, @std, @median}); 
    [x_meanVal_vigor,x_maxVal_vigor, x_stdVal_vigor, x_medVal_vigor] = grpstats(choiceforcerate_table.abs_dUtil_normed,Y_vigor,{@mean, @max, @std, @median}); 

    figure(520+nn)
    hold on

    tblnew = table();
    tblnew.abs_dUtil_normed = domain';
    tblnew.subj = 0*ones(size(tblnew.abs_dUtil_normed));
    tblnew.subj = categorical(tblnew.subj);
    [test,test_CI] = predict(mdl_delta_util,tblnew);
    plot(domain,test,'Color',[150,100,150]./255,'LineWidth',3);
    x_CI = [domain,domain(end:-1:1)];
    y_CI = [test_CI(:,1)',test_CI(end:-1:1,2)'];
    p = fill(x_CI,y_CI,'red');p.FaceColor = [150,100,150]./255;p.EdgeColor = 'None';p.FaceAlpha = 0.5;

    errorbar(x_meanVal_vigor,meanVal_vigor,stdVal_vigor/sqrt(nsubjs),stdVal_vigor/sqrt(nsubjs),...
        x_stdVal_vigor./sqrt(nsubjs),x_stdVal_vigor./sqrt(nsubjs),'k','Linewidth',2,'LineStyle','None')%SE
%     plot(domain,mdl_delta_util.Coefficients(1,2).Estimate + mdl_delta_util.Coefficients(2,2).Estimate*domain,'Color',[150,100,150]./255,'LineWidth',3)

    xlabel('Normalized Difference in Subjective Utility')
    ylabel('Choice Vigor - Rate of Force, Normalized to Max (%MVC/sec)')
    title('(B) Choice Vigor vs. Difference in Utility')

    AICs = [AICs,mdl_delta_util.ModelCriterion.AIC];
    BICs = [BICs,mdl_delta_util.ModelCriterion.BIC];
    Ftests = [Ftests,coefTest(mdl_delta_util)];
    MDLS{nn} = mdl_delta_util;
end

%% Print AICS
AICs
%% Print Model Parameters

%% Dont Plot Losing Models
[~,ind_winner] = min(AICs);

if 1 
    for nn = 1:length(modelstocompare)
        if nn ~= ind_winner
            close(figure(520+nn))
        end
    end
end
