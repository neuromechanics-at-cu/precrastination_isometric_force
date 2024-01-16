%% Example Choice Behavior
% Figure 1 in Manuscript
% DELIBERATION TIME AND MAX CURSOR VEL

if ~exist('example_trial')
    load('../../Data/Processed Data/processed_data.mat','example_trial')
end

figure(101)
subplot(2,1,1);hold on;title('Cursor Position (px)')
% plot(Mchoice{2}.framedata(739).time,Mchoice{2}.framedata(739).x_choice_px)
% plot(Mchoice{2}.framedata(268).time,Mchoice{2}.framedata(268).x_choice_px)
% for github
plot(example_trial.time,example_trial.x_choice_px)
% legend('Low Diff in Utility','High Diff in Utility')
ylim([440,840])
xlim([0,0.85])
subplot(2,1,2);hold on;title('Cursor Velocity (px/sec)')
% plot(Mchoice{2}.framedata(739).time(1:end-1),Mchoice{2}.framedata(739).dx_px_choice_smooth)
% plot(Mchoice{2}.framedata(268).time(1:end-1),Mchoice{2}.framedata(268).dx_px_choice_smooth)
% for git hub
plot(example_trial.time(1:end-1),example_trial.dx_px_choice_smooth)
% legend('Low Diff in Utility','High Diff in Utility')
ylim([-1000 5000])
xlim([0,0.85])