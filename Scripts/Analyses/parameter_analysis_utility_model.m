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


%% Report Parameter Values Across All Models
param_table = nan(4,6,2);
for ii = 1:4 % Parameter
    for jj = 1:6 % Model
        param_table(ii,jj,:) = bootci(10000,@mean,squeeze(sol(ii,jj,:)));
    end
end


%% Compare Parameters of Winning Model for CoVariance

figure
combos = nchoosek(1:4,2);

var_str = {'Temp. Sens. Effort', 'Temp. Sens. Task', 'log(Temperature)', 'Weighting'};

sol_temp = sol;
sol_temp(3,:,:) = log(sol_temp(3,:,:));

for ii = 1:size(combos,1)
%     plot(squeeze(sol(combos(ii,1),4,:)),squeeze(sol(combos(ii,2),4,:)),'x')
%     tbl = array2table([squeeze(sol(combos(ii,1),4,:))';squeeze(sol(combos(ii,2),4,:))']');
    tbl = array2table([squeeze(sol_temp(combos(ii,1),4,:))';squeeze(sol_temp(combos(ii,2),4,:))']');
    mdl = fitlm(tbl,'Var2~Var1')
    pvalues(ii) = mdl.Coefficients.pValue(2);
    subplot(3,3,(combos(ii,1)-1)*3 + combos(ii,2)-1)
    plot(mdl)
    ylabel(var_str{combos(ii,2)})
    xlabel(var_str{combos(ii,1)})
    title([var_str{combos(ii,1)},' vs. ',var_str{combos(ii,2)}])
end

% And if we do this, then nothing is significant...
bootci(10000,@nanmean,squeeze(sol(1,bestmodel,:)))
bootci(10000,@nanmean,squeeze(sol(2,bestmodel,:)))
bootci(10000,@nanmean,squeeze(sol(4,bestmodel,:)))

