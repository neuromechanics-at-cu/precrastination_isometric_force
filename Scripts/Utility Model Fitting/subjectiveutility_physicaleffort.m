function subj_u = subjectiveutility_physicaleffort(gamma,reward,totalduration,rampuptime2base,rampuptime2top,...
                                        rampdowntime2base,rampdowntime2end,baselevel,toplevel,...
                                        durat,onset,discounttype)

% Time step
% dt = 0.1;% sec %CMH NEW 2023
dt = 0.01;% sec %CMH NEW 2023
% dt = 0.001;% sec
% dt = 0.00001;% sec


% % In case onset is not a multiple of dt
% onset = onset - mod(onset,dt);


% forcetimeintegral(0.01,17,0.5,0.5,0.5,0.5,0.1,0.5,1,2.5);
baseline = 225;

timevec = 0:dt:totalduration;

% Effort
[~,f_vsTime] = forcetimeintegral(dt,totalduration,rampuptime2base,rampuptime2top,...
                                            rampdowntime2base,rampdowntime2end,baselevel,toplevel,...
                                            durat,onset);

% Reward ... reward should be 1000
rewardvec = ones(size(timevec))*reward/length(timevec);
                                        
                                        
% Net Utility
if length(rewardvec) ~= length(f_vsTime)
%     obj_uvec = rewardvec - f_vsTime(1:end-1);
    obj_uvec = rewardvec - f_vsTime(1:end-1)/baseline;
else
%     obj_uvec = rewardvec - f_vsTime;
    obj_uvec = rewardvec - f_vsTime/baseline;
end                    
                        
% Subjective Utilty
if discounttype == 0 % Exponential
    subj_uvec = obj_uvec.*exp(-1*gamma*timevec);
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
