function [m] = processData_Familiarization(path,subjs,mvc,el_ang)


if nargin < 4
    disp('ERROR: Missing input arguments.')
    return
end

numsubjs = length(subjs);

for s = 1:numsubjs
    subj2analyze = subjs{s};    
    datadir = [path,'/',subj2analyze];

    %% Import Robot Data for given subject
    [T]=robotdataread_Effort_Familiarization(subj2analyze,datadir);

    T.accuracy = zeros(T.trials,1);
    for ii = 1:T.trials
        %% Calculate Accuracy
        T.framedata(ii).y_profile_px = nan(size(T.framedata(ii).adc2));
        for jj = 1:length(T.framedata(ii).adc2)
            T.framedata(ii).y_profile_px(jj) = adc2ypx(T.framedata(ii).adc2(jj),T.framedata(ii).adc4(jj),...
                                                T.framedata(ii).adc2(1),T.framedata(ii).adc4(1),...
                                                el_ang(s),mvc(s));
        end
        startind = find(T.framedata(ii).statenumber == 5,1,'first');
        endind = find(T.framedata(ii).statenumber == 5,1,'last');
        for jj = startind:endind
            T.accuracy(ii) = T.accuracy(ii) + abs(T.framedata(ii).y_profile_px(jj) - profileFromParams(T.framedata(ii).cursorx(jj),T.trialduration(ii),T.rampuptime2base(ii),...
                                                                                            T.rampuptime2top(ii),T.rampdowntime2base(ii),T.rampdowntime2end(ii),...
                                                                                            T.baselevel(ii),T.toplevel(ii),T.topduration(ii),T.toponset(ii)));
        end
        T.accuracy(ii) = T.accuracy(ii)/(endind - startind + 1);
        
        %% Calculate Peak Rate of Force "Vigor"
        % Calculate Force Trace
        T.framedata(ii).y_force = -cosd(270 - el_ang(s))*T.framedata(ii).ft_xworld + sind(270 - el_ang(s))*T.framedata(ii).ft_zworld;
        %TODO: T.x_force
%         coeff = ones(1,34)/34;
%         y_smooth = filtfilt(coeff,1,T.y_profile_px)';
        fs = 200;%Hz Sampling freq
        fc = 10;%Hz cutoff freq
        Wn = fc/(fs/2);
        [b,a] = butter(4,Wn,'low');
        y_smooth = filtfilt(b,a,T.framedata(ii).y_profile_px);
        dy_smooth = diff(y_smooth)./diff(T.framedata(ii).time);
        % Constrain to area of force generation for hill or valley
        windowstart = startind+T.toponset(ii)*200 - 100;
        windowend = startind+(T.toponset(ii)+T.rampuptime2top(ii)+T.topduration(ii)+T.rampdowntime2base(ii))*200 + 100;
        T.peakRateOfForceUphill_px(ii,1) = -1*min(dy_smooth(windowstart:windowend));
        T.peakRateOfForceDownhill_px(ii,1) = -1*max(dy_smooth(windowstart:windowend));
%         dF_smooth = filtfilt(coeff,1,T.y_force)';
        F_smooth = filtfilt(b,a,T.framedata(ii).y_force);
        dF_smooth = diff(F_smooth)./diff(T.framedata(ii).time);
        T.peakRateOfForceUphill_N(ii,1) = -1*min(dF_smooth(windowstart:windowend));
        T.peakRateOfForceDownhill_N(ii,1) = -1*max(dF_smooth(windowstart:windowend));
    end
    
    %% Record any failed trials
    T.failedTrials = [];
    for jj = 1:T.trials
        if find(T.framedata(jj).avstatenumber == 9)
            T.failedTrials = [T.failedtrials,jj];
        end        
    end
    
    %% Accuracy Summary Calcs
    % Calculate average accuracy by trial type over the familiarization
    % period. Would be better to just look at last trial for each type
    % though, due to learning. Min trial type is 0, thus indexing +1.
    for ii = min(T.trialtype):max(T.trialtype)
        ind = find(T.trialtype == ii,1);
        [T.forcetimeintegral_trialtype(ii+1),~] = forcetimeintegral(0.001,T.trialduration(ind),T.rampuptime2base(ind),T.rampuptime2top(ind),T.rampdowntime2base(ind),...
                                                    T.rampdowntime2end(ind),T.baselevel(ind),T.toplevel(ind),T.topduration(ind),T.toponset(ind));
    end
    
    for ii = min(T.trialtype):max(T.trialtype)
        T.mean_accuracytot(ii+1) = mean(T.accuracy(find(T.trialtype == ii)));
    end
    % Lump by Duration
    durats = unique(T.topduration);
    for ii = 1:length(durats)
        T.mean_accuracyByDuration(ii) = mean(T.accuracy(find(T.topduration == durats(ii))));
    end
    % Lump by Onset Time
    onsets = unique(T.toponset);
    for ii = 1:length(onsets)
        T.mean_accuracyByOnset(ii) = mean(T.accuracy(find(T.toponset == onsets(ii))));
    end
    % Lump by Off Time
    T.offtime = T.toponset + T.rampuptime2top + T.topduration + T.rampdowntime2base;
    offtime = unique(T.offtime);
    for ii = 1:length(offtime)
        T.mean_accuracyByOfftime(ii) = mean(T.accuracy(find(T.offtime == offtime(ii))));
    end
    % Lump by distance to center
    T.offcenter = min(abs(T.trialduration/2 - T.toponset),abs(T.trialduration/2 - T.offtime));
    T.offcenter2 = abs((T.toponset + (T.rampuptime2top + T.topduration + T.rampdowntime2base)/2) - T.trialduration/2);
    
    offcenter = unique(T.offcenter);
    for ii = 1:length(offcenter)
        T.mean_accuracyByOffcenter(ii) = mean(T.accuracy(find(T.offcenter == offcenter(ii))));
    end
    offcenter2 = unique(T.offcenter2);
    for ii = 1:length(offcenter2)
        T.mean_accuracyByOffcenter2(ii) = mean(T.accuracy(find(T.offcenter2 == offcenter2(ii))));
    end
    %Just last trial (of 3 familiarization trials)
    for ii = min(T.trialtype):max(T.trialtype)
        T.mean_accuracytot_last1(ii+1) = mean(T.accuracy(find(T.trialtype(1:3*(max(T.trialtype)+1)) == ii,1,'last')));
    end
    % Lump by Duration
    for ii = 1:length(durats)
        temp = 3:3:3*(max(T.trialtype)+1);
        ind = find(T.topduration(3:3:3*(max(T.trialtype)+1)) == durats(ii));
        T.mean_accuracyByDuration_last1(ii) = mean(T.accuracy(temp(ind)));
    end
    % Lump by Onset Time
    for ii = 1:length(onsets)
        temp = 3:3:3*(max(T.trialtype)+1);
        ind = find(T.toponset(3:3:3*(max(T.trialtype)+1)) == onsets(ii));
        T.mean_accuracyByOnset_last1(ii) = mean(T.accuracy(temp(ind)));
    end
%     % Just last trials when randomly presented
%     %Just last couple trials
%     for ii = min(T.trialtype):max(T.trialtype)
%         T.mean_accuracytot_lastrand(ii+1) = mean(T.accuracy(find(T.trialtype == ii,1,'last')));
%     end
%     % Lump by Duration
%     for ii = 1:length(durats)
%         T.mean_accuracyByDuration_lastrand(ii) = mean(T.accuracy(find(T.topduration(3*(max(T.trialtype)+1)+1:end) == durats(ii)) + 3*(max(T.trialtype)+1)));
%     end
%     % Lump by Onset Time
%     for ii = 1:length(onsets)
%         T.mean_accuracyByOnset_lastrand(ii) = mean(T.accuracy(find(T.toponset(3*(max(T.trialtype)+1)+1:end) == onsets(ii)) + 3*(max(T.trialtype)+1)));
%     end

    % Just last trials when randomly presented
    %Just last couple trials
    for ii = min(T.trialtype):max(T.trialtype)
        T.mean_accuracytot_lastrand(ii+1) = mean(T.accuracy(find(T.trialtype == ii,2,'last')));
    end
    % Lump by Duration
    for ii = 1:length(durats)
        T.mean_accuracyByDuration_lastrand(ii) = mean(T.accuracy(find(T.topduration(3*(max(T.trialtype)+1)+1:end) == durats(ii)) + 3*(max(T.trialtype)+1)));
    end
    % Lump by Onset Time
    for ii = 1:length(onsets)
        T.mean_accuracyByOnset_lastrand(ii) = mean(T.accuracy(find(T.toponset(3*(max(T.trialtype)+1)+1:end) == onsets(ii)) + 3*(max(T.trialtype)+1)));
    end

m{s} = T;
end
