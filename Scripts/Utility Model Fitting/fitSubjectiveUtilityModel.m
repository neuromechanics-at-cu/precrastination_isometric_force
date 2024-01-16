function cost = fitSubjectiveUtilityModel(x,Mchoice,reward,utilitycombotype,costtype,taskmodeltype,discounttype)
% Fit Just Raw Data

% Task Model Type
% 0 - Fdot squared
% 1 - Carry until hill/valley complete (binary)
% 2 - Binary hill or valley

% Set/Initialize Params
beta = x(3);
% p_predicted = nan(length(Mchoice.trials));
p_predicted = nan(1,Mchoice.trials);


% Choose parameters based on model type ("Utility Combo Type")
if utilitycombotype == 1  % No Temporal Sensitivity, Force Only
    rel_weight = 1;
    gamma_physical = 0;
    gamma_task = 0;
elseif utilitycombotype == 2 % Temporal Sensitivity, Force Only
    rel_weight = 1;
    gamma_physical = x(1);
%     gamma_task = x(1);
    gamma_task = 0;
elseif utilitycombotype == 3% No temporal Sensitivity, Combined
    rel_weight = x(4);
    gamma_physical = 0;    
    gamma_task = 0;
elseif utilitycombotype == 4 % Temporal Sensitivity, Combined
    rel_weight = x(4);
    gamma_physical = x(1);
    gamma_task = x(1);
elseif utilitycombotype == 5 % Temporal Sensitivity, Combined, Separate Gammas
    rel_weight = x(4);
    gamma_physical = x(1);
    gamma_task = x(2);
elseif utilitycombotype == 6 % Temporal Sensitivity Fdot only, TN Force 
    rel_weight = x(4);
    gamma_physical = 0;
    gamma_task = x(2);
elseif utilitycombotype == 7 % Time in State Only
    rel_weight = -1;
    gamma_physical = 0;
    gamma_task = x(2);
    taskmodeltype = 8;
elseif utilitycombotype == 8 % Temporal Sensitivity F + Time in State 
    rel_weight = x(4);
    gamma_physical = x(1);
    gamma_task = x(2);
    taskmodeltype = 8;
else
    error('Unspecified type.')
end

%% Find estimated probability for each choice
% based on parameters and subjective utility
% Note this could be done more efficiently if done by choice type, rather
% than all choices

for ii = 1:Mchoice.trials
    if rel_weight == 1
        subj_u_L = subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice.Ltrialduration(ii),Mchoice.Lrampuptime2base(ii),Mchoice.Lrampuptime2top(ii),...
                            Mchoice.Lrampdowntime2base(ii),Mchoice.Lrampdowntime2end(ii),Mchoice.Lbaselevel(ii),Mchoice.Ltoplevel(ii),...
                            Mchoice.Ltopduration(ii),Mchoice.Ltoponset(ii),discounttype);    
        subj_u_R = subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice.Rtrialduration(ii),Mchoice.Rrampuptime2base(ii),Mchoice.Rrampuptime2top(ii),...
                            Mchoice.Rrampdowntime2base(ii),Mchoice.Rrampdowntime2end(ii),Mchoice.Rbaselevel(ii),Mchoice.Rtoplevel(ii),...
                            Mchoice.Rtopduration(ii),Mchoice.Rtoponset(ii),discounttype); 
    elseif rel_weight == -1 % Hack - Task Only
        subj_u_L = subjectiveutility_taskeffort(gamma_task,reward,Mchoice.Ltrialduration(ii),Mchoice.Lrampuptime2base(ii),Mchoice.Lrampuptime2top(ii),...
                                Mchoice.Lrampdowntime2base(ii),Mchoice.Lrampdowntime2end(ii),Mchoice.Lbaselevel(ii),Mchoice.Ltoplevel(ii),...
                                Mchoice.Ltopduration(ii),Mchoice.Ltoponset(ii),taskmodeltype,discounttype);
        subj_u_R = subjectiveutility_task(gamma_task,reward,Mchoice.Rtrialduration(ii),Mchoice.Rrampuptime2base(ii),Mchoice.Rrampuptime2top(ii),...
                                Mchoice.Rrampdowntime2base(ii),Mchoice.Rrampdowntime2end(ii),Mchoice.Rbaselevel(ii),Mchoice.Rtoplevel(ii),...
                                Mchoice.Rtopduration(ii),Mchoice.Rtoponset(ii),taskmodeltype,discounttype); 

    else
        subj_u_L = rel_weight*subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice.Ltrialduration(ii),Mchoice.Lrampuptime2base(ii),Mchoice.Lrampuptime2top(ii),...
                                Mchoice.Lrampdowntime2base(ii),Mchoice.Lrampdowntime2end(ii),Mchoice.Lbaselevel(ii),Mchoice.Ltoplevel(ii),...
                                Mchoice.Ltopduration(ii),Mchoice.Ltoponset(ii),discounttype) + ...
           (1-rel_weight)*subjectiveutility_task(gamma_task,reward,Mchoice.Ltrialduration(ii),Mchoice.Lrampuptime2base(ii),Mchoice.Lrampuptime2top(ii),...
                                Mchoice.Lrampdowntime2base(ii),Mchoice.Lrampdowntime2end(ii),Mchoice.Lbaselevel(ii),Mchoice.Ltoplevel(ii),...
                                Mchoice.Ltopduration(ii),Mchoice.Ltoponset(ii),taskmodeltype,discounttype);
        subj_u_R = rel_weight*subjectiveutility_physicaleffort(gamma_physical,reward,Mchoice.Rtrialduration(ii),Mchoice.Rrampuptime2base(ii),Mchoice.Rrampuptime2top(ii),...
                                Mchoice.Rrampdowntime2base(ii),Mchoice.Rrampdowntime2end(ii),Mchoice.Rbaselevel(ii),Mchoice.Rtoplevel(ii),...
                                Mchoice.Rtopduration(ii),Mchoice.Rtoponset(ii),discounttype) + ...
           (1-rel_weight)*subjectiveutility_task(gamma_task,reward,Mchoice.Rtrialduration(ii),Mchoice.Rrampuptime2base(ii),Mchoice.Rrampuptime2top(ii),...
                                Mchoice.Rrampdowntime2base(ii),Mchoice.Rrampdowntime2end(ii),Mchoice.Rbaselevel(ii),Mchoice.Rtoplevel(ii),...
                                Mchoice.Rtopduration(ii),Mchoice.Rtoponset(ii),taskmodeltype,discounttype); 
    end
    % Probability of Choosing Left Choice
    p_predicted(ii) = ( 1 + exp(-1*beta*(subj_u_L - subj_u_R)) ).^-1;                    
end

% Prevent log(0) NaN Error
p_predicted(p_predicted==0) = eps;
p_predicted(p_predicted==1) = 1-eps;
% p_predicted(p_predicted==0) = 2*eps;
% p_predicted(p_predicted==1) = 1+2*eps;

% Determine Quality of Fit / Cost to be minimized
% RMS Error
choseLeft = -1*Mchoice.choice'+2;
goodTrials = ~isnan(choseLeft);
if costtype == 1
    cost = sqrt(sum((choseLeft(goodTrials) - p_predicted(goodTrials)).^2)/length(p_predicted(goodTrials)));
% Negative Log Likelihood
elseif costtype == 2 
    cost = -1*sum(log(p_predicted(goodTrials)).*choseLeft(goodTrials) + log(1-p_predicted(goodTrials)).*(1-choseLeft(goodTrials))); 
%     p_pred_no_zero = p_predicted;
%     p_pred_no_zero(p_pred_no_zero==0) = eps;
%     cost_term_1 = log(p_pred_no_zero(goodTrials)).*choseLeft(goodTrials);
% 
%     p_pred_no_one = p_predicted;
%     p_pred_no_one(p_pred_no_one==1) = 1-eps;
%     cost_term_2 = log(1-p_pred_no_one(goodTrials)).*(1-choseLeft(goodTrials));
% 
%     cost = -1*sum(cost_term_1+cost_term_2); 

else
    error('Unspecified cost type.')
end