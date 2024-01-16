%% Parameter Analysis for Winning Model
% Table 3 in Manuscript

if ~exist('sol')
    load('../Utility Model Fitting/utility_model_fit_results.mat','sol')
end

%%%
bestmodel = 4; %from pxp in compare_utility_models.m
%%%

% Physical Effort Time Sensitivity
bootci(10000,@mean,squeeze(sol(1,bestmodel,:)))
% Task Cost Time Sensitivity
bootci(10000,@mean,squeeze(sol(2,bestmodel,:)))
% Temperature
bootci(10000,@mean,squeeze(sol(3,bestmodel,:)))
% Proportion Task vs. Physical Effort
bootci(10000,@mean,squeeze(sol(4,bestmodel,:)))
