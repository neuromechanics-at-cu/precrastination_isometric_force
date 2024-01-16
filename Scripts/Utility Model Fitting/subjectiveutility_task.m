function subj_u = subjectiveutility_task(gamma,reward,totalduration,rampuptime2base,rampuptime2top,...
                                        rampdowntime2base,rampdowntime2end,baselevel,toplevel,...
                                        durat,onset,modeltype,discounttype)

% Time step
% dt = 0.1;% sec % CMH NEW 2023
dt = 0.01;% sec % CMH NEW 2023
% dt = 0.001;% sec
% dt = 0.00001;% sec

% [baseline,~] = dforce_dt_timeintegral(dt,17,0.5,0.5,0.5,0.5,0.1,0.5,...
%                                             1,2.5);
baseline = 68.04; %shouldnt it be 69?

% % In case onset is not a multiple of dt
% onset = onset - mod(onset,dt);

timevec = 0:dt:totalduration;
taskload_vs_time = [];

% Task Carry Effort 
if modeltype == 0 % dF_dt
    [~,taskload_vs_time] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,...
                                            rampdowntime2base,rampdowntime2end,baselevel,toplevel,...
                                            durat,onset);
elseif modeltype == 1 % carry load until complete
    offset_time = onset + durat + rampuptime2top + rampdowntime2base;
    taskload_vs_time = [ones(1,length(0:dt:offset_time)-1),zeros(1,length(offset_time:dt:totalduration))];
elseif modeltype == 2 % load during hill or valley
    offset_time = onset + durat + rampuptime2top + rampdowntime2base;
    taskload_vs_time = [zeros(1,length(0:dt:onset)-1),...
                  ones(1,length(onset:dt:offset_time)-1),...
                  zeros(1,length(offset_time:dt:totalduration))];
elseif modeltype == 3 % Fdot, but binary
    taskload_vs_time = [ones(1,length(0:dt:rampuptime2base)-1),...
                        zeros(1,length(rampuptime2base:dt:onset)-1),...
                        ones(1,length(onset:dt:onset+rampuptime2top)-1),...
                        zeros(1,length(onset+rampuptime2top:dt:onset+rampuptime2top+durat)-1),...
                        ones(1,length(onset+rampuptime2top+durat:dt:onset+rampuptime2top+durat+rampdowntime2base)-1),...
                        zeros(1,length(onset+rampuptime2top+durat+rampdowntime2base:dt:totalduration-rampdowntime2end)-1),...
                        ones(1,length(totalduration-rampdowntime2end:dt:totalduration))];
elseif modeltype == 4 % Ramp up and down and full hill/valley
    taskload_vs_time = [ones(1,length(0:dt:rampuptime2base)-1),...
                        zeros(1,length(rampuptime2base:dt:onset)-1),...
                        ones(1,length(onset:dt:onset+rampuptime2top+durat+rampdowntime2base)-1),...
                        zeros(1,length(onset+rampuptime2top+durat+rampdowntime2base:dt:totalduration-rampdowntime2end)-1),...
                        ones(1,length(totalduration-rampdowntime2end:dt:totalduration))];
elseif modeltype == 5 % Carry until Start of Hill or Valley
    taskload_vs_time = [ones(1,length(0:dt:onset)-1),zeros(1,length(onset:dt:totalduration))];
elseif modeltype == 6 % Binary, "Blips" At beginning and end of just Hill/Valley
    taskload_vs_time = [zeros(1,length(0:dt:onset)-1),...
                        ones(1,length(onset:dt:onset+rampuptime2top)-1),...
                        zeros(1,length(onset+rampuptime2top:dt:onset+rampuptime2top+durat)-1),...
                        ones(1,length(onset+rampuptime2top+durat:dt:onset+rampuptime2top+durat+rampdowntime2base)-1),...
                        zeros(1,length(onset+rampuptime2top+durat+rampdowntime2base:dt:totalduration))];
elseif modeltype == 7 % Fdot, but only for force generation (upslopes)
    [~,taskload_vs_time] = dforce_dt_timeintegral(dt,totalduration,rampuptime2base,rampuptime2top,...
                                            rampdowntime2base,rampdowntime2end,baselevel,toplevel,...
                                            durat,onset);
    taskload_vs_time = taskload_vs_time.*(taskload_vs_time>0);

elseif modeltype == 8 % Time in State / Constant Force
    [~,taskload_vs_time] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durat,onset);
    timevec = [zeros(1,length(0:dt:(rampuptime2base-dt))),... % rampup time, not penalized
                0:dt:onset-rampuptime2base-dt,... % Hill 1
                zeros(1,length(0:dt:(rampuptime2top-dt))),... % rampdown time, not penalized
                0:dt:durat-dt,... % valley
                zeros(1,length(0:dt:(rampdowntime2base-dt))),... % rampuptime, not penalized
                0:dt:(totalduration-onset-rampuptime2top-rampdowntime2base-durat-rampdowntime2end-dt),... % hill 2
                zeros(1,length(0:dt:(rampdowntime2end)))];
end

% Reward ... reward should be 1000
rewardvec = ones(size(timevec))*reward/length(timevec);                                     
                                        
% Net Utility
if modeltype == 8 
    obj_uvec = rewardvec - taskload_vs_time; % Dont square for the case.
else
    if length(rewardvec) ~= length(taskload_vs_time)
%         obj_uvec = rewardvec - taskload_vs_time(1:end-1).^2;
        obj_uvec = rewardvec - (taskload_vs_time(1:end-1).^2)/baseline;
    else
        obj_uvec = rewardvec - (taskload_vs_time.^2)/baseline;
    end
end

% Subjective Utilty
if discounttype == 0 % Exponential
%     subj_uvec = obj_uvec.*exp(-1*gamma*timevec);
elseif discounttype == 1 % Hyperbolic
%     subj_uvec = obj_uvec.*(1./(sign(gamma)+gamma*timevec) - (sign(gamma)-1)/2);
    subj_uvec = obj_uvec.*1/(1+gamma*timevec);
% elseif discounttype = 2 % Logistic Growth
%     subj_uvec = 1./(1+exp(-1*gamma*timevec)); 
elseif discounttype == 2 % Exponential, Equally Weighted 2023
    if gamma < 0
        subj_uvec = obj_uvec.*exp(-1*gamma*timevec);
    else
        subj_uvec = obj_uvec.*exp(gamma*(timevec(end)-timevec));
    end   
end
subj_u = sum(subj_uvec);
