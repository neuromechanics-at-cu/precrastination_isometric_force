%% Load In Data
load('../../Data/Processed Data/processed_data.mat')

%% Do Subjects Fatigue? Is MVC correlated to test setup/elbow position?
% Investigate MVC within subject, from beginning and end, session to
% session

mvc_A = [mvc1A;mvc2A;mvc3A];
mvc_B = [mvc1B;mvc2B;mvc3B];
mvc_1_se = std([mvc1A,mvc1B])/sqrt(length(subjs));
mvc_2_se = std([mvc2A,mvc2B])/sqrt(length(subjs));
mvc_3_se = nanstd([mvc3A,mvc3B])/sqrt(length(subjs));
mvc_1_mean = mean([mvc1A,mvc1B]);
mvc_2_mean = mean([mvc2A,mvc2B]);
mvc_3_mean = nanmean([mvc3A,mvc3B]);

% %Plot mvc beginning and end
% figure
% hold on
% plot([1,2,3],mvc_A,'Color',[0.8,0.8,0.8])
% plot([1,2,3],mvc_B,'Color',[0.8,0.8,0.8])
% p = plot([1,2,3],[mvc_1_mean,mvc_2_mean,mvc_3_mean]./mvc_1_mean,'LineWidth',3,'Color',[0,0,0]);
% e = errorbar([1,2,3],[mvc_1_mean,mvc_2_mean,mvc_3_mean],[mvc_1_se,mvc_2_se,mvc_3_se]*1.96);
% set(e,'LineWidth',3,'Color',[0,0,0])
% axis([0.8 3.2 0 5])
% ylabel('MVC')
% xticks([1 2 3])
% xticklabels({'Session Start','End of Fam/Start of Choice','Session End'})


mvc_A_norm = [mvc1A;mvc2A;mvc3A]./mvc1A;
mvc_B_norm = [mvc1B;mvc2B;mvc3B]./mvc1B;
mvc_1_se = nanstd([mvc_A_norm(1,:),mvc_B_norm(1,:)])/sqrt(length(subjs));
mvc_2_se = nanstd([mvc_A_norm(2,:),mvc_B_norm(2,:)])/sqrt(length(subjs));
mvc_3_se = nanstd([mvc_A_norm(3,:),mvc_B_norm(3,:)])/sqrt(length(subjs));
mvc_1_mean = nanmean([mvc_A_norm(1,:),mvc_B_norm(1,:)]);
mvc_2_mean = nanmean([mvc_A_norm(2,:),mvc_B_norm(2,:)]);
mvc_3_mean = nanmean([mvc_A_norm(3,:),mvc_B_norm(3,:)]);
figure(601)
hold on
plot([1,2,3],mvc_A_norm,'Color',[0.8,0.8,0.8])
plot([1,2,3],mvc_B_norm,'Color',[0.8,0.8,0.8])
p = plot([1,2,3],[mvc_1_mean,mvc_2_mean,mvc_3_mean],'LineWidth',3,'Color',[0,0,0]);
e = errorbar([1,2,3],[mvc_1_mean,mvc_2_mean,mvc_3_mean],[mvc_1_se,mvc_2_se,mvc_3_se]*1.96);
set(e,'LineWidth',3,'Color',[0,0,0])
axis([0.8 3.2 0.5 1.5])
ylabel('MVC')
xticks([1 2 3])
xticklabels({'Session Start','End of Fam/Start of Choice','Session End'})


%% Subjective Fatigue

fatigue_A = [fatigue1A;fatigue2A;fatigue3A];
fatigue_B = [fatigue1B;fatigue2B;fatigue3B];
fatigue_1_se = std([fatigue1A,fatigue1B])/sqrt(length(subjs));
fatigue_2_se = std([fatigue2A,fatigue2B])/sqrt(length(subjs));
fatigue_3_se = nanstd([fatigue3A,fatigue3B])/sqrt(length(subjs));
fatigue_1_mean = mean([fatigue1A,fatigue1B]);
fatigue_2_mean = mean([fatigue2A,fatigue2B]);
fatigue_3_mean = nanmean([fatigue3A,fatigue3B]);

%Plot subjective fatigue beginning and end
figure(602)
hold on
plot([1,2,3],fatigue_A,'Color',[0.8,0.8,0.8])
plot([1,2,3],fatigue_B,'Color',[0.8,0.8,0.8])
p = plot([1,2,3],[fatigue_1_mean,fatigue_2_mean,fatigue_3_mean],'LineWidth',3,'Color',[0,0,0]);
e = errorbar([1,2,3],[fatigue_1_mean,fatigue_2_mean,fatigue_3_mean],[fatigue_1_se,fatigue_2_se,fatigue_3_se]*1.96);
set(e,'LineWidth',3,'Color',[0,0,0])
axis([0.8 3.2 0 5])
ylabel('Subjective Fatigue (Likert, 1-5)')
xticks([1 2 3])
xticklabels({'Session Start','End of Fam/Start of Choice','Session End'})



% % Make a Big Table
stages = {'Beginning','Post-Fam/Pre-Choice','End'};
session_types = {'Hill','Valley'};
fatigue_comb = [fatigue_A,fatigue_B];
mvc_comb = [mvc_A,mvc_B];
norm_mvc_comb = [mvc_A./mvc1A,mvc_B./mvc1B];
newmat = nan(length(subjs)*length(session_types)*length(stages),5);
for ii = 1:length(subjs)
    for jj = 1:length(session_types)
        for kk = 1:length(stages)
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),1) = ii;
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),2) = jj;
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),3) = kk;
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),4) = fatigue_comb(kk,ii+(jj-1)*length(subjs));
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),5) = mvc_comb(kk,ii+(jj-1)*length(subjs));
            newmat(kk+(jj-1)*length(stages)+(ii-1)*length(session_types)*length(stages),6) = norm_mvc_comb(kk,ii+(jj-1)*length(subjs));
        end
    end
end

bigtbl = array2table(newmat,'VariableNames',{'subjnum','session_type','stage','fatigue_rating','mvc','normed_mvc'});

bigtbl.subjnum = categorical(bigtbl.subjnum);
bigtbl.session_type = categorical(bigtbl.session_type);
bigtbl.stage = ordinal(bigtbl.stage);
% bigtbl.fatigue_rating = ordinal(bigtbl.fatigue_rating);


%% Linear Model of MVC (Objective Fatigue) 
fitlme(bigtbl,'mvc ~ stage + session_type + (1|subjnum)')


%% Lienar Model of Fatigue Rating 
fitlme(bigtbl,'fatigue_rating ~ stage + session_type + (1|subjnum)')


