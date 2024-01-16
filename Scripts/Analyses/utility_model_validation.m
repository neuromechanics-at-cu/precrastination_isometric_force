
if ~exist('AIC')
    run('compare_utility_models.m')
end

addpath('../Utility Model Fitting/')
taskmodeltype = 0; 
discounttype = 2;

%% Calculate Subjective Utility Based on Model Fits

modelstocompare = 4;
reward = 0;


include_outlier_subjs = 1;
subj_outliers = [];

for nn = 1:length(modelstocompare)
    selectedmodel = modelstocompare(nn)*ones(length(subjs),1);

    % Based on Best Fit Model Structure for Majority of Subjects
    % selectedmodel = modelstocompare(find(bmsout.pxp==max(bmsout.pxp)))*ones(length(subjs),1);
    % OR Based on Best Fit Model for that Subject
    % selectedmodel = bestmodel;
    
    sol_reduced = squeeze(sol); % should be (model, subj, param)
    
%     subj_outliers = isoutlier(AIC(:,modelstocompare(nn)));
    subj_outliers = isoutlier(AIC(:,modelstocompare(nn)),'percentile',[5 95]);

    if include_outlier_subjs == 0
        nsubjs = sum(~subj_outliers);
    else
        nsubjs = length(subjs);
    end

    choice_table_prob = [];
    % For Each Subject
    for jj = 1:length(subjs)
        if include_outlier_subjs == 1 || (include_outlier_subjs == 0 && subj_outliers(jj) == 0)
        
            gamma_physical = sol_reduced(1,selectedmodel(jj),jj);
            gamma_task = sol_reduced(2,selectedmodel(jj),jj);
            beta = sol_reduced(3,selectedmodel(jj),jj);
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
            Mchoice{jj}.delta_utility_minimaxed = (Mchoice{jj}.delta_utility - min(abs(Mchoice{jj}.delta_utility)))/(max(abs(Mchoice{jj}.delta_utility))-min(abs(Mchoice{jj}.delta_utility)));
            %% Prepare Data Table
            for ii = 1:Mchoice{jj}.trials
                if ~isnan(Mchoice{jj}.choice(ii))
%                     choice_table_prob = [choice_table_prob;
%                                         jj,Mchoice{jj}.delta_utility(ii),Mchoice{jj}.delta_utility_normed(ii) ,abs(Mchoice{jj}.abs_delta_utility_minimaxed(ii)),Mchoice{jj}.delta_utility_minimaxed(ii),...
%                                            Mchoice{jj}.p_predicted(ii),-1*Mchoice{jj}.choice(ii)+2];
                    choice_table_prob = [choice_table_prob;
                                        jj,ii>455,Mchoice{jj}.delta_utility(ii),Mchoice{jj}.delta_utility_normed(ii),abs(Mchoice{jj}.abs_delta_utility_minimaxed(ii)),Mchoice{jj}.delta_utility_minimaxed(ii),...
                                           Mchoice{jj}.p_predicted(ii),-1*Mchoice{jj}.choice(ii)+2,...
                                           Mchoice{jj}.subj_u_L(ii),Mchoice{jj}.Ltoponset(ii),Mchoice{jj}.Ltopduration(ii),...
                                           Mchoice{jj}.subj_u_R(ii),Mchoice{jj}.Rtoponset(ii),Mchoice{jj}.Rtopduration(ii)];
                end
            end
        end
    end    

    %% Label Data Table
    choice_table_prob = array2table(choice_table_prob,'VariableNames',{'subj','session_type','subj_dUtil','dUtil_normed','abs_dUtil_normed','dUtil_minimaxed',...
                                                                       'p_choice','chose_left',...
                                                                       'subjUtil_L','onset_L','duration_L',...
                                                                       'subjUtil_R','onset_R','duration_R'});
    choice_table_prob.subj = categorical(choice_table_prob.subj);


    %% Bin and Average For Plotting...
    bin_midpts = -1:.2:1;
    bins_for_model = -1:0.1:1;
    [Y,E] = discretize(choice_table_prob.dUtil_minimaxed,bin_midpts);
    [Y_model,E] = discretize(choice_table_prob.dUtil_minimaxed,bins_for_model);
    
    [meanVal, maxVal, stdVal, medVal] = grpstats(choice_table_prob.p_choice,Y_model,{@mean, @max, @std, @median}); 
    [x_meanVal,x_maxVal, x_stdVal, x_medVal] = grpstats(choice_table_prob.dUtil_minimaxed,Y_model,{@mean, @max, @std, @median}); 

    figure(401); hold on;
%     errorbar(x_meanVal,meanVal,stdVal/sqrt(nsubjs),stdVal/sqrt(nsubjs),'k','Linewidth',2)%SE
    % MODEL 
    x_CI = [bins_for_model,bins_for_model(end-1:-1:2)];
    y_CI = [meanVal'+1.96*stdVal'/sqrt(nsubjs),meanVal(end:-1:1)'-1.96*stdVal(end:-1:1)'/sqrt(nsubjs)];
    p = fill(x_CI,y_CI,'red');p.FaceColor = [150,100,150]./255;p.EdgeColor = 'None';p.FaceAlpha = 0.2;
    plot(x_meanVal,meanVal,'Color',[150,100,150]./255,'Linewidth',2)%SE
%     errorbar(x_meanVal,meanVal,stdVal/sqrt(nsubjs),stdVal/sqrt(nsubjs),'k','Linewidth',2)%SE

    hold on
    [Ncounts,E2] = histcounts(choice_table_prob.dUtil_minimaxed,bin_midpts);
    [data_meanVal, data_maxVal, data_stdVal, data_medVal] = grpstats(choice_table_prob.chose_left,Y,{@mean, @max, @std, @median,}); 
    [data_x_meanVal, data_x_maxVal, data_x_stdVal, data_x_medVal] = grpstats(choice_table_prob.dUtil_minimaxed,Y,{@mean, @max, @std, @median,}); 
%     errorbar(data_x_meanVal,data_meanVal,data_x_stdVal/sqrt(nsubjs),data_stdVal/sqrt(nsubjs),'Color',[150,100,150]./255,'Linewidth',2)
    errorbar(data_x_meanVal,data_meanVal,data_stdVal/sqrt(nsubjs),data_stdVal/sqrt(nsubjs),'Color','k','Linewidth',2)

end
    xlim([-1 1])
    ylim([0 1])
    xlabel('Normalized Difference in Subjective Utility')
    ylabel('p(Choosing Left)')
    title('(A) Average Model Prediction')

%% Generate Model Predictions Based on Average Parameter Values


%% Plot Example Subject....
% Subject 10
% Task Procrastinator
% Force PreCrastinator

% Phys \ Task | Pre | Pro |
%        Pre  |  7  |  9  |
%        Pro  |  5  |  4  |

% The most frequent is Physical Pre, Task Pro
figure(402)
hold on

subjtable = choice_table_prob(choice_table_prob.subj == categorical(6),:);

unique_combos = unique(subjtable(:,{'subj_dUtil','p_choice'}));
sorted_unique_combos = sortrows(unique_combos,'subj_dUtil');
% plot(sorted_unique_combos.dUtil_minimaxed,sorted_unique_combos.p_choice)
plot(sorted_unique_combos.subj_dUtil,sorted_unique_combos.p_choice,':','Color',[150,100,150]./255,'LineWidth',2)

unique_combos_data = unique(subjtable(:,{'subj_dUtil','chose_left'}));

tblstats = grpstats(subjtable,"subj_dUtil","mean","DataVars","chose_left");
plot(tblstats.subj_dUtil,tblstats.mean_chose_left,'kx')

%
% bin_midpts = -1.25:.5:1.25;
bin_midpts = -1.5:.5:1.5;
[Y,E] = discretize(subjtable.subj_dUtil,bin_midpts);
[data_meanVal, data_maxVal, data_stdVal, data_medVal] = grpstats(subjtable.chose_left,Y,{@mean, @max, @std, @median}); 
[data_x_meanVal, data_x_maxVal, data_x_stdVal, data_x_medVal] = grpstats(subjtable.subj_dUtil,Y,{@mean, @max, @std, @median}); 
% errorbar(data_x_meanVal,data_meanVal,data_x_stdVal/sqrt(nsubjs),data_stdVal/sqrt(nsubjs),'Color','k','Linewidth',2)
errorbar(data_x_meanVal,data_meanVal,data_stdVal/sqrt(nsubjs),data_stdVal/sqrt(nsubjs),'Color','k','Linewidth',2)
xlim([-1.25 1.25])
ylim([0 1])
title('(B) Subject 6 - Example Preference Curve')
xlabel('Difference in Subjective Utility')
ylabel('p(Choosing Left)')



% Do it for Subj 6
ii = 6;
% ii = 18;
% ii = 2;

subjtable = choice_table_prob(choice_table_prob.subj == categorical(ii),:);

% Hack to Make Up for my onesided thing...
LeftProfiles  = unique(subjtable(:,{'subjUtil_L','duration_L','onset_L','session_type'})); 
LeftProfiles = renamevars(LeftProfiles,{'subjUtil_L','duration_L','onset_L'},{'subjUtil','duration','onset'});
RightProfiles = unique(subjtable(:,{'subjUtil_R','duration_R','onset_R','session_type'}));
RightProfiles = renamevars(RightProfiles,{'subjUtil_R','duration_R','onset_R'},{'subjUtil','duration','onset'});
%         Right_and_Left = join(LeftProfiles,RightProfiles);
[~,ind] = setdiff(RightProfiles.subjUtil,LeftProfiles.subjUtil);
Right_and_Left = [LeftProfiles;RightProfiles(ind,:)];

unique_util = unique(Right_and_Left);
unique_util = sortrows(unique_util,'subjUtil','descend');
unique_util.rank = [1:28]';

% figure;
% plot(1:28,unique(Mchoice{2}.subj_u_L),'x')

figure(403); hold on;
session_color = {'b','r'};
duration_symbol = {'x','+','o'};
durs = 1:2:5;
sess = [0,1];
for kk = 1:length(sess)
    for mm = 1:length(durs)% duration
        criteria = (unique_util.session_type == sess(kk)) & (unique_util.duration == durs(mm));
        plot(unique_util.rank(criteria),unique_util.subjUtil(criteria),'Marker',duration_symbol(mm),'Color',session_color{kk})
    end
end
        legend({'Hills','Valleys'})
ylabel("Subjective Utility")
xlabel('Rank from Worst to Best')
title(['(C) - Subj: ',num2str(ii),' Ranked Profiles'])