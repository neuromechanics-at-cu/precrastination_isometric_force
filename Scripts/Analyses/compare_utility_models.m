%% Compare Models of Subjective Utility
% Found in TABLE 2 of Manuscript

addpath('../Supporting Functions/')

load('../Utility Model Fitting/utility_model_fit_results.mat')
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end

AIC = nan(length(modeltypes(modelstocompare)),length(subjs));
BIC = nan(length(modeltypes(modelstocompare)),length(subjs));
Likelihood = -1*fval;
for ii = 1:length(modeltypes(modelstocompare))
    for jj = 1:length(subjs)
        AIC(ii,jj) = 2*freeparams(ii) - 2*Likelihood(ii,jj);
        BIC(ii,jj) = Likelihood(ii,jj) - freeparams(ii)*log(sum(~isnan(Mchoice{jj}.choice)))/2;
    end
end
AIC = -1*AIC'; % Multiply this by -1. SPM uses a negative BIC, where max is better.
[bmsout.alpha,bmsout.exp_r,bmsout.xp,bmsout.pxp,bmsout.bor] = spm_BMS(AIC,[],1);

figure(21)
subplot(2,2,1)
h(1)=bar(bmsout.alpha);
set(gca,'XTickLabel',modelnames(modelstocompare))
title('alpha')
subplot(2,2,2)
h(2)=bar(bmsout.exp_r);
set(gca,'XTickLabel',modelnames(modelstocompare))
title('exp_r')
subplot(2,2,3)
h(3)=bar(bmsout.xp);
set(gca,'XTickLabel',modelnames(modelstocompare))
title('xp')
subplot(2,2,4)
h(4)=bar(bmsout.pxp);
set(gca,'XTickLabel',modelnames(modelstocompare))
title('pxp')

figure(22)
h(5)=bar(bmsout.pxp);
set(gca,'XTickLabel',modelnames(modelstocompare))
title('pxp')
disp(modelnames(modelstocompare)')
disp('Exceedance Probability')
disp(bmsout.xp)
disp('Protected Exceedance Probability')
disp(bmsout.pxp)
disp('BOR')
disp(bmsout.bor)

lme_mean=mean(AIC)
lme_se=std(AIC)./sqrt(length(subjs))
disp([lme_mean-1.96*lme_se; lme_mean+1.96*lme_se ])


for jj = 1:length(subjs)
    [~,ind(jj)] = max(AIC(jj,:));
end
