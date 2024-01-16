%% Load In Variables
% subjs
% sol
% Mchoice
if ~exist('Mchoice')
    load('../../Data/Processed Data/processed_data.mat','Mchoice')
end
if ~exist('subjs')
    load('../../Data/Processed Data/processed_data.mat','subjs')
end
if ~exist('sol')
    load('../Utility Model Fitting/utility_model_fit_results.mat','sol')
end
% Other Vars
taskmodeltype = 0; 
discounttype = 2;
reward = 0;

include_outlier_trials = 0; 
include_outlier_subjs = 1;
modelstocompare = 1:4;

addpath('../Utility Model Fitting/')
addpath('../Supporting Functions/')


%% Deliberation Time vs Utility
deliberation_time_analysis

%% Choice Vigor vs Utility
choice_vigor_analysis

%% Deliberation Time vs. Choice Vigor
choice_vigor_vs_delib_time