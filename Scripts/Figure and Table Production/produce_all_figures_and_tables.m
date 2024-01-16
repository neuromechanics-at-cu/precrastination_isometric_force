%% Figure 1 - Example Choice Metrics
example_choice_behaviors
%% Figure 2 - Example Subjective Utility of Profiles
example_utility_costs
%% Table 1 - Linear Model of Preference Results
run('../Analyses/linearmodel_preference.m')
%% Figure 3 - Linear Model of Preference Plots
run('../Analyses/plot_linearmodel_preference.m')
%% Table 2 - Subjective Utility Model Comparison
run('../Analyses/compare_utility_models.m')
%% Figure 4 - Subjective Utility Model Summary and Example Subject 
run('../Analyses/utility_model_validation.m')
%% Table 3 - Bootstrapping Subjective Utility Model Parameters
run('../Analyses/parameter_analysis_utility_model.m')
%% Table 4 - Choice Behavior Linear Model Results, and 
%% Figure 5 - Choice Behavior Model Plots
run('../Analyses/choice_behavior_analysis.m')
%% Figure 6 - Fatigue Models and Plots
run('../Analyses/fatigue_analysis.m')
%% Other Analyses
% Accuracy
run('../Analyses/accuracy_analysis.m')
% Failed and Discarded Trials
run('../Analyses/failed_and_discarded_trials.m')