%% Plot Demonstration of Subjective Utility Model Costs
% Figure 2 in Manuscript
% Compare Early/Late Hills and Valleys, Time Sensitivity, and Cumulative Costs

addpath('../Supporting Functions/')

%% Set Discounting Factor
gamma = 0.1;

%% Initialize Hill Parameters
dt = 0.001;
% Force Level
baselevel = 0.1;
toplevel = 0.5;
% Times
rampuptime2base = 0.5;% (s)
rampuptime2top = 0.5;% (s)
rampdowntime2base = 0.5;% (s)
rampdowntime2end = 0.5;
totalduration = 17;% (s)
onsets = [2.5:2:12.5;
2.5,4.5,8.5,10.5,NaN,NaN;
2.5:2:8.5,NaN,NaN];
durats = [1,3,5]; % (s), Choose #2, 3 sec
time = 0:dt:totalduration;

%% EARLY HILL

%%% Force 
[~,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,1));
Uforce_TD = f_vsTime.*(exp(-gamma*time));

figure(201)
subplot(2,2,1)
hold on;
plot(time,f_vsTime);
xlim([0 totalduration])
ylim([0 0.55])
title('(A) Physical Effort Cost - Early/Late Short Hill')

%%% Rate of Force
[~,df_dt_vsTime] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,1));

%%% Rate of Force Squared
Uforcedotsq = (df_dt_vsTime.^2).*exp(gamma*time)/exp(gamma*time(end));
figure(201)
subplot(2,2,3)
hold on;
plot(time,df_dt_vsTime.^2);
xlim([0 totalduration])
ylim([0 1])
title('(C) Task Cost - Early/Late Short Hill')

%% LATE HILL

%%% Force 
[~,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,6));
Uforce_TD = f_vsTime.*(exp(-gamma*time));

figure(201)
subplot(2,2,1)
hold on;
plot(time,f_vsTime);
xlim([0 totalduration])
ylim([0 0.55])

%%% Rate of Force
[~,df_dt_vsTime] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,6));

%%% Rate of Force Squared
Uforcedotsq = (df_dt_vsTime.^2).*exp(gamma*time)/exp(gamma*time(end));
figure(201)
subplot(2,2,3)
hold on;
plot(time,df_dt_vsTime.^2);



%% Initialize Valley Parameters
% Force Level
baselevel = 0.2515;
toplevel = 0.05;
% Times
rampuptime2base = 1;% (s)
rampuptime2top = 0.5;% (s)
rampdowntime2base = 0.5;% (s)
rampdowntime2end = 1;
totalduration = 17;% (s)
onsets = [2.5:2:12.5;
2.5,4.5,8.5,10.5,NaN,NaN;
2.5:2:8.5,NaN,NaN];
durats = [1,3,5]; % (s), Choose #2, 3 sec
time = 0:dt:totalduration;

%% EARLY VALLEY
%%% Force 
[~,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,1));
Uforce_TD = f_vsTime.*(exp(-gamma*time));


figure(201)
subplot(2,2,2)
hold on;
plot(time,f_vsTime);
xlim([0 totalduration])
ylim([0 0.55])
title('(B) Physical Effort Cost - Early/Late Short Valley')

%%% Rate of Force
[~,df_dt_vsTime] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,1));

%%% Rate of Force Squared
Uforcedotsq = (df_dt_vsTime.^2).*exp(gamma*time)/exp(gamma*time(end));
figure(201)
subplot(2,2,4)
hold on;
plot(time,df_dt_vsTime.^2);
plot(time,Uforcedotsq)
xlim([0 totalduration])
ylim([0 1])
title('(D) Task Cost - Early/Late Short Valley')


% Plot Combo
alpha = 0.5;
figure(203);
subplot(1,2,1)
plot(time,alpha*Uforcedotsq/68.04 + (1-alpha)*Uforce_TD/225)
hold on;
plot(time,0.5*Uforcedotsq/68.04)
plot(time,0.5*Uforce_TD/225) 
cumcost(1) = sum(0.5*Uforcedotsq/68.04) + sum(0.5*Uforce_TD/225);
xlim([0 17])
title('(G) Early Valley')


%% LATE VALLEY
%%% Force
[fti_exact,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,6));
Uforce_TD = f_vsTime.*(exp(-gamma*time));


figure(201)
subplot(2,2,2)
hold on;
plot(time,f_vsTime);
xlim([0 totalduration])
ylim([0 0.55])

%%% Rate of Force
[~,df_dt_vsTime] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(1),onsets(1,6));

%%% Rate of Force Squared
Uforcedotsq = (df_dt_vsTime.^2).*exp(gamma*time)/exp(gamma*time(end));
figure(201)
subplot(2,2,4)
hold on;
plot(time,df_dt_vsTime.^2);
xlim([0 totalduration])
ylim([0 1])

alpha = 0.5;
figure(203)
subplot(1,2,2)
plot(time,alpha*Uforcedotsq/68.04 + (1-alpha)*Uforce_TD/225)
hold on;
plot(time,0.5*Uforcedotsq/68.04)
plot(time,0.5*Uforce_TD/225) 
cumcost(2) = sum(0.5*Uforcedotsq/68.04) + sum(0.5*Uforce_TD/225);
xlim([0 17])
title('(H) Late Valley')


%% Time Sensitivity
figure(202)
hold on
plot(time,exp(-gamma*time))
plot(time,exp(gamma*time)/exp(gamma*time(end)))
ylim([0 1])
xlim([0 17])
xlabel('Time (s)')
ylabel('Scaling Factor')
title('(E) Time Sensitivity')

%% Bar Plot
figure(204)
bar(categorical({'Early Valley','Late Valley'}),cumcost)
ylabel('Cumulative Cost (a.u.)')
title('(F)')