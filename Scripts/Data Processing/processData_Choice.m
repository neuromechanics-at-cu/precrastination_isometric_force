function [m] = processData_Choice(path,subjs,mvc,el_ang)

if nargin < 4
    disp('ERROR: Missing input arguments.')
    return
end

numsubjs = length(subjs);

for s = 1:numsubjs
    subj2analyze = subjs{s};    
    datadir = [path,'/',subj2analyze];

    %% Import Robot Data for given subject
    [T]=robotdataread_Effort_Choice(subj2analyze,datadir);
    trialdatanames=fieldnames(T);
    framedatanames=fieldnames(T.framedata);
    configdatanames=fieldnames(T.config);
    
    T.choice = zeros(T.trials,1);
    T.choicetime = zeros(T.trials,1);
    T.accuracy = zeros(T.trials,1);
    %% Find choices, reaction time
    for ii = 1:T.trials
        % Check if Bad Trial
%         if T.framedata(ii).time(end) > 0.1
        if T.framedata(ii).time(end) > 0.2
            if T.realization(ii)
                ind = find(T.framedata(ii).statenumber == 3,1,'last');
%                 if T.framedata(ii).cursorx(ind) > 700
                if T.framedata(ii).cursorx(ind) > 670
                    T.choice(ii) = 2; % Right
                elseif T.framedata(ii).cursorx(ind) < 610
%                 elseif T.framedata(ii).cursorx(ind) < 580
                    T.choice(ii) = 1; % Left
                else
                    T.choice(ii) = NaN; % Just in case...
                end
%             elseif T.framedata(ii).cursorx(end) > 700
            elseif T.framedata(ii).cursorx(end) > 670
                T.choice(ii) = 2;
%             elseif T.framedata(ii).cursorx(end) < 580
            elseif T.framedata(ii).cursorx(end) < 610
                T.choice(ii) = 1;
            else 
                T.choice(ii) = NaN;
            end
            % Choice Time (sec)
            T.choicetime(ii) = (find(T.framedata(ii).statenumber == 3,1,'last') - find(T.framedata(ii).statenumber == 3,1,'first'))/200;

            % Force Time Integral
            [T.Lforcetimeintegral(ii),~] = forcetimeintegral(0.001,T.Ltrialduration(ii),T.Lrampuptime2base(ii),T.Lrampuptime2top(ii),T.Lrampdowntime2base(ii),...
                                                            T.Lrampdowntime2end(ii),T.Lbaselevel(ii),T.Ltoplevel(ii),T.Ltopduration(ii),T.Ltoponset(ii));
            [T.Rforcetimeintegral(ii),~] = forcetimeintegral(0.001,T.Rtrialduration(ii),T.Rrampuptime2base(ii),T.Rrampuptime2top(ii),T.Rrampdowntime2base(ii),T.Rrampdowntime2end(ii),...
                                                            T.Rbaselevel(ii),T.Rtoplevel(ii),T.Rtopduration(ii),T.Rtoponset(ii));
            %% Force Trace For All Trials
            T.framedata(ii).y_force = -cosd(270 - el_ang(s))*T.framedata(ii).ft_xworld + sind(270 - el_ang(s))*T.framedata(ii).ft_zworld;
            T.framedata(ii).x_force = -sind(270 - el_ang(s))*T.framedata(ii).ft_xworld + cosd(270 - el_ang(s))*T.framedata(ii).ft_zworld;
            T.framedata(ii).tot_force = (T.framedata(ii).x_force.^2 + T.framedata(ii).y_force.^2).^0.5;
            T.framedata(ii).dx_force_dt_smooth = diff(T.framedata(ii).x_force)./diff(T.framedata(ii).time);

            %% X Pixel Trace, Rate of Force of Choice (Vigor Estimation)
            T.framedata(ii).x_choice_px = nan(size(T.framedata(ii).adc2));
            x_px_pre1 = 1280/2;
            x_px_pre2 = 1280/2;
            for jj = 1:length(T.framedata(ii).adc2)

                T.framedata(ii).x_choice_px(jj) = adc2xpx(T.framedata(ii).adc2(jj),T.framedata(ii).adc4(jj),...
                                                    T.framedata(ii).adc2(1),T.framedata(ii).adc4(1),...
                                                    el_ang(s),mvc(s),x_px_pre1,x_px_pre2);
                x_px_pre2 = x_px_pre1;
                x_px_pre1 = T.framedata(ii).x_choice_px(jj);
            end
            % Smooth Pixels, Find Rate
            fs = 200;%Hz Sampling freq
            fc = 10;%Hz cutoff freq
            Wn = fc/(fs/2);
            [b,a] = butter(4,Wn,'low');
            x_px_smooth = filtfilt(b,a,T.framedata(ii).x_choice_px);
            T.framedata(ii).dx_px_choice_smooth = diff(x_px_smooth)./diff(T.framedata(ii).time);
            % Find Max X Px Rate
            startind = find(T.framedata(ii).statenumber == 3,1,'first');
            endind = min([find(T.framedata(ii).statenumber == 3,1,'last'),length(T.framedata(ii).dx_px_choice_smooth)]);
            if ~isempty(startind) && endind > 0 
                T.max_dxpx_dt(ii) = max(abs(T.framedata(ii).dx_px_choice_smooth(startind:endind)));
            else
                T.max_dxpx_dt(ii) = NaN;
            end

            % Find Deliberation Time (sec) based on velocity threshold of 200 px/s
            enddelib = find(abs(T.framedata(ii).dx_px_choice_smooth) > 200,1,'first');
            if ~isempty(enddelib)
                T.delibtime(ii) = T.framedata(ii).time(enddelib);
            else
                T.delibtime(ii) = NaN;
            end

            %% Accuracy, if realized
            if T.realization(ii)
                T.framedata(ii).y_profile_px = nan(size(T.framedata(ii).adc2));
                for jj = 1:length(T.framedata(ii).adc2)
                    T.framedata(ii).y_profile_px(jj) = adc2ypx(T.framedata(ii).adc2(jj),T.framedata(ii).adc4(jj),...
                                                        T.framedata(ii).adc2(1),T.framedata(ii).adc4(1),...
                                                        el_ang(s),mvc(s));
                end
                % Calculate Force Trace
%                 T.framedata(ii).y_force = -cosd(270 - el_ang(s))*T.framedata(ii).ft_xworld + sind(270 - el_ang(s))*T.framedata(ii).ft_zworld;
                startind = find(T.framedata(ii).statenumber == 6,1,'first');
                endind = find(T.framedata(ii).statenumber == 7,1,'first');
                y_smooth = filtfilt(b,a,T.framedata(ii).y_profile_px);
                dy_smooth = diff(y_smooth)./diff(T.framedata(ii).time);


                y_smooth_subset = filtfilt(b,a,T.framedata(ii).y_profile_px(startind:endind));
                mvc_subset = (1024- y_smooth_subset - 100)/824;%convert px to %mvc
                time_subset = T.framedata(ii).time(startind:endind) - T.framedata(ii).time(startind);
                % Find Peak Rate of Force
                rate_mvc_subset = diff(mvc_subset)./diff(time_subset);


                % Calculate Peak Rate of Force "Vigor"
                if T.choice(ii) == 1 % left
%                     coeff = ones(1,34)/34;
%                     y_smooth = filtfilt(coeff,1,T.y_profile_px)';
%                     fs = 200;%Hz Sampling freq
%                     fc = 10;%Hz cutoff freq
%                     Wn = fc/(fs/2);
%                     [b,a] = butter(4,Wn,'low');
%                     y_smooth = filtfilt(b,a,T.framedata(ii).y_profile_px);
%                     dy_smooth = diff(y_smooth)./diff(T.framedata(ii).time);
                    % Constrain to area of force generation for hill or valley
                    windowstart = startind+T.Ltoponset(ii)*200 - 100;
                    windowend = startind+(T.Ltoponset(ii)+T.Lrampuptime2top(ii)+T.Ltopduration(ii)+T.Lrampdowntime2base(ii))*200 + 100;
                    T.peakRateOfForceUphill_px(ii) = -1*min(dy_smooth(windowstart:windowend));
                    T.peakRateOfForceDownhill_px(ii) = -1*max(dy_smooth(windowstart:windowend));
%                     dF_smooth = filtfilt(coeff,1,T.y_force)';
                    F_smooth = filtfilt(b,a,T.framedata(ii).y_force);
                    dF_smooth = diff(F_smooth)./diff(T.framedata(ii).time);
                    T.peakRateOfForceUphill_N(ii) = -1*min(dF_smooth(windowstart:windowend));
                    T.peakRateOfForceDownhill_N(ii) = -1*max(dF_smooth(windowstart:windowend));
                    % New Rev 1,
                    % Base 1 Indices
                    [~,base1_startind] = min(abs(time_subset - T.Lrampuptime2base(ii)));
                    [~,base1_endind] = min(abs(time_subset - T.Ltoponset(ii)));

                    % Hill / Valley Indices
                    [~,hillvalley_startind] = min(abs(time_subset - (T.Ltoponset(ii) + T.Lrampuptime2top(ii))));
                    [~,hillvalley_endind] = min(abs(time_subset - (T.Ltoponset(ii) + T.Lrampuptime2top(ii) + T.Ltopduration(ii))));

                    % Base 2 Indicies
                    [~,base2_startind] = min(abs(time_subset - (T.Ltoponset(ii) + T.Lrampuptime2top(ii) + T.Ltopduration(ii) + T.Lrampdowntime2base(ii))));
                    [~,base2_endind] = min(abs(time_subset - (T.Ltrialduration(ii) - T.Lrampdowntime2end(ii))));

                    % Calculate Mean Forces 
                    % Base 1
                    T.avg_force_base1(ii) = mean(mvc_subset(base1_startind:base1_endind));
                    T.std_force_base1(ii) = std(mvc_subset(base1_startind:base1_endind));
                    % Hill/Valley
                    T.avg_force_hillvalley(ii) = mean(mvc_subset(hillvalley_startind:hillvalley_endind));
                    T.std_force_hillvalley(ii) = std(mvc_subset(hillvalley_startind:hillvalley_endind));
                    % Base 2
                    T.avg_force_base2(ii) = mean(mvc_subset(base2_startind:base2_endind));
                    T.std_force_base2(ii) = std(mvc_subset(base2_startind:base2_endind));

                    % Find Overshoot and Peak Rate, Ramp  +/- 1 second
                    if T.Ltoplevel(ii) == 0.5 % Hill
                        % Ramp 1 (uphill)
                        T.overshoot_ramp1(ii) = max(mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        T.peakrate_mvc_ramp1(ii) = max(rate_mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        
                        % Ramp 2 (downhill)
                        T.overshoot_ramp2(ii) = min(mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                        T.peakrate_mvc_ramp2(ii) = min(rate_mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                    else % Valley
                        % Ramp 1 (downhill)
                        T.overshoot_ramp1(ii) = min(mvc_subset((base1_endind-200):(hillvalley_startind+200)));            
                        T.peakrate_mvc_ramp1(ii) = min(rate_mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        % Ramp 2 (uphill)
                        T.overshoot_ramp2(ii) = max(mvc_subset((hillvalley_startind-200):(base2_startind+200)));                
                        T.peakrate_mvc_ramp2(ii) = max(rate_mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                    end


                elseif T.choice(ii) == 2 % right
%                     coeff = ones(1,34)/34;
%                     y_smooth = filtfilt(coeff,1,T.y_profile_px)';
%                     fs = 200;%Hz Sampling freq
%                     fc = 10;%Hz cutoff freq
%                     Wn = fc/(fs/2);
%                     [b,a] = butter(4,Wn,'low');
%                     y_smooth = filtfilt(b,a,T.framedata(ii).y_profile_px);
%                     dy_smooth = diff(y_smooth)./diff(T.framedata(ii).time);
                    % Constrain to area of force generation for hill or valley
                    windowstart = startind+T.Rtoponset(ii)*200 - 100;
                    windowend = startind+(T.Rtoponset(ii)+T.Rrampuptime2top(ii)+T.Rtopduration(ii)+T.Rrampdowntime2base(ii))*200 + 100;
                    T.peakRateOfForceUphill_px(ii) = -1*min(dy_smooth(windowstart:windowend));
                    T.peakRateOfForceDownhill_px(ii) = -1*max(dy_smooth(windowstart:windowend));
%                     dF_smooth = filtfilt(coeff,1,T.y_force)';
                    F_smooth = filtfilt(b,a,T.framedata(ii).y_force);
                    dF_smooth = diff(F_smooth)./diff(T.framedata(ii).time);
                    T.peakRateOfForceUphill_N(ii) = -1*min(dF_smooth(windowstart:windowend));
                    T.peakRateOfForceDownhill_N(ii) = -1*max(dF_smooth(windowstart:windowend));

                    % New Rev 1,
                    % Base 1 Indices
                    [~,base1_startind] = min(abs(time_subset - T.Rrampuptime2base(ii)));
                    [~,base1_endind] = min(abs(time_subset - T.Rtoponset(ii)));

                    % Hill / Valley Indices
                    [~,hillvalley_startind] = min(abs(time_subset - (T.Rtoponset(ii) + T.Rrampuptime2top(ii))));
                    [~,hillvalley_endind] = min(abs(time_subset - (T.Rtoponset(ii) + T.Rrampuptime2top(ii) + T.Rtopduration(ii))));

                    % Base 2 Indicies
                    [~,base2_startind] = min(abs(time_subset - (T.Rtoponset(ii) + T.Rrampuptime2top(ii) + T.Rtopduration(ii) + T.Rrampdowntime2base(ii))));
                    [~,base2_endind] = min(abs(time_subset - (T.Rtrialduration(ii) - T.Rrampdowntime2end(ii))));
                    % Calculate Mean Forces 
                    % Base 1
                    T.avg_force_base1(ii) = mean(mvc_subset(base1_startind:base1_endind));
                    T.std_force_base1(ii) = std(mvc_subset(base1_startind:base1_endind));
                    % Hill/Valley
                    T.avg_force_hillvalley(ii) = mean(mvc_subset(hillvalley_startind:hillvalley_endind));
                    T.std_force_hillvalley(ii) = std(mvc_subset(hillvalley_startind:hillvalley_endind));
                    % Base 2
                    T.avg_force_base2(ii) = mean(mvc_subset(base2_startind:base2_endind));
                    T.std_force_base2(ii) = std(mvc_subset(base2_startind:base2_endind));
                    % Find Overshoot and Peak Rate, Ramp  +/- 1 second
                    if T.Rtoplevel(ii) == 0.5 % Hill
                        % Ramp 1 (uphill)
                        T.overshoot_ramp1(ii) = max(mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        T.peakrate_mvc_ramp1(ii) = max(rate_mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        
                        % Ramp 2 (downhill)
                        T.overshoot_ramp2(ii) = min(mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                        T.peakrate_mvc_ramp2(ii) = min(rate_mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                    else % Valley
                        % Ramp 1 (downhill)
                        T.overshoot_ramp1(ii) = min(mvc_subset((base1_endind-200):(hillvalley_startind+200)));            
                        T.peakrate_mvc_ramp1(ii) = min(rate_mvc_subset((base1_endind-200):(hillvalley_startind+200)));
                        % Ramp 2 (uphill)
                        T.overshoot_ramp2(ii) = max(mvc_subset((hillvalley_startind-200):(base2_startind+200)));                
                        T.peakrate_mvc_ramp2(ii) = max(rate_mvc_subset((hillvalley_startind-200):(base2_startind+200)));
                    end


                else
                    T.peakRateOfForceUphill_px(ii) = NaN;
                    T.peakRateOfForceDownhill_px(ii) = NaN;
                    T.peakRateOfForceUphill_N(ii) = NaN;
                    T.peakRateOfForceDownhill_N(ii) = NaN;
                    T.avg_force_base1(ii) = NaN;
                    T.std_force_base1(ii) = NaN;
                    T.avg_force_hillvalley(ii) = NaN;
                    T.std_force_hillvalley(ii) = NaN;
                    T.avg_force_base2(ii) = NaN;
                    T.std_force_base2(ii) = NaN;
                    T.overshoot_ramp1(ii) = NaN;
                    T.peakrate_mvc_ramp1(ii) = NaN;
                    T.overshoot_ramp2(ii) = NaN;
                    T.peakrate_mvc_ramp2(ii) = NaN;

                end
                

                for jj = startind:endind
                    if T.choice(ii) == 1 % Left
%                         T.accuracy(ii) = T.accuracy(ii) + abs(recreateY(T.framedata(ii).adc2(jj),T.framedata(ii).adc4(jj),...
%                                                                 T.framedata(ii).adc2(1),T.framedata(ii).adc4(1),...
%                                                                 234,mvc(s)) - profileFromParams(T.framedata(ii).cursorx(jj),T.Ltrialduration(ii),T.Lrampuptime(ii),...
%                                                                                                 T.Lrampdowntime(ii),T.Lbaselevel(ii),T.Ltoplevel(ii),T.Ltopduration(ii),T.Ltoponset(ii)));
                        T.accuracy(ii) = T.accuracy(ii) + abs(T.framedata(ii).y_profile_px(jj) - profileFromParams(T.framedata(ii).cursorx(jj),T.Ltrialduration(ii),T.Lrampuptime2base(ii),...
                                                                                                T.Lrampuptime2top(ii),T.Lrampdowntime2base(ii),T.Lrampdowntime2end(ii),...
                                                                                                T.Lbaselevel(ii),T.Ltoplevel(ii),T.Ltopduration(ii),T.Ltoponset(ii)));
                    elseif T.choice(ii) == 2 %Right
                        T.accuracy(ii) = T.accuracy(ii) + abs(T.framedata(ii).y_profile_px(jj) - profileFromParams(T.framedata(ii).cursorx(jj),T.Rtrialduration(ii),T.Rrampuptime2base(ii),...
                                                                                                T.Rrampuptime2top(ii),T.Rrampdowntime2base(ii),T.Rrampdowntime2end(ii),...
                                                                                                T.Rbaselevel(ii),T.Rtoplevel(ii),T.Rtopduration(ii),T.Rtoponset(ii)));
                    else 
                        T.accuracy(ii) = NaN;
                    end
                end
                T.accuracy(ii) = T.accuracy(ii)/(endind - startind + 1);
            else
                T.accuracy(ii) = NaN;
                T.peakRateOfForceUphill_px(ii) = NaN;
                T.peakRateOfForceDownhill_px(ii) = NaN;
                T.peakRateOfForceUphill_N(ii) = NaN;
                T.peakRateOfForceDownhill_N(ii) = NaN;
                % Add More Performance Measures
                T.avg_force_base1(ii) = NaN;
                T.std_force_base1(ii) = NaN;
                T.avg_force_hillvalley(ii) = NaN;
                T.std_force_hillvalley(ii) = NaN;
                T.avg_force_base2(ii) = NaN;
                T.std_force_base2(ii) = NaN;
    
                T.overshoot_ramp1(ii) = NaN;
                T.peakrate_mvc_ramp1(ii) = NaN;
                T.overshoot_ramp2(ii) = NaN;
                T.peakrate_mvc_ramp2(ii) = NaN;
            end
        else
            T.choice(ii) = NaN;
            T.choicetime(ii) = NaN;
            T.accuracy(ii) = NaN;
            T.peakRateOfForceUphill_px(ii) = NaN;
            T.peakRateOfForceDownhill_px(ii) = NaN;
            T.peakRateOfForceUphill_N(ii) = NaN;
            T.peakRateOfForceDownhill_N(ii) = NaN;
            % Add More Performance Measures
            T.avg_force_base1(ii) = NaN;
            T.std_force_base1(ii) = NaN;
            T.avg_force_hillvalley(ii) = NaN;
            T.std_force_hillvalley(ii) = NaN;
            T.avg_force_base2(ii) = NaN;
            T.std_force_base2(ii) = NaN;

            T.overshoot_ramp1(ii) = NaN;
            T.peakrate_mvc_ramp1(ii) = NaN;
            T.overshoot_ramp2(ii) = NaN;
            T.peakrate_mvc_ramp2(ii) = NaN;
        end
    end

    %% Normalize, Remove Outliers for Delib and Choice Time
    T.delibtime_nooutliers = T.delibtime;
    T.delibtime_nooutliers(isoutlier(T.delibtime,'percentiles',[1 99])) = NaN;
    T.mean_delibtime_nooutliers = nanmean(T.delibtime_nooutliers);
    T.delibtime_nooutliers_normed = T.delibtime_nooutliers - T.mean_delibtime_nooutliers;
    % CMH NEW
    T.delibtime_nooutliers_zcored = (T.delibtime_nooutliers - T.mean_delibtime_nooutliers)./nanstd(T.delibtime_nooutliers);

    T.mean_delibtime = nanmean(T.delibtime);
    T.delibtime_normed = T.delibtime - T.mean_delibtime;

    %% Record any failed trials
    T.failedTrials = [];
    for jj = 1:T.trials
        if find(T.framedata(jj).avstatenumber == 9)
            T.failedTrials = [T.failedtrials,jj];
        end        
    end
        
    %% Look at Choices
    %TODO: Revisit these
    % Removing for now, until I'm certain we need this. 2021.

    numTimes = max(T.choicetypecount); %5
    numChoices = T.trials/numTimes; %455/5 = 91
    numTypes = max([T.Ltrialtype',T.Rtrialtype'])+1;% [14]
    combos = nchoosek([0:numTypes-1],2);

    % 3 duration conditions, 5 start conditions
    durats = unique(T.Ltopduration);
    onsets = unique(T.Ltoponset);
    T.choiceMatrix = zeros(length(durats),length(onsets));
    T.choiceArray = zeros(1,numTypes);
    nDurats = zeros(size(durats));
    nOnsets = zeros(size(onsets));
    T.choiceTypeMatrix = zeros(numChoices,2);
    T.combotype = nan(T.trials,1); % Initialize (NEW)

    % Accumulate choices in matrix of conditions...
    for ii = 1:T.trials
        if T.choice(ii) == 1
            row = find(T.Ltopduration(ii) == durats);
            col = find(T.Ltoponset(ii) == onsets);
            T.choiceArray(T.Ltrialtype(ii)+1) = T.choiceArray(T.Ltrialtype(ii)+1) + 1;
            T.choiceMatrix(row,col) = T.choiceMatrix(row,col) + 1;
        elseif T.choice(ii) == 2
            row = find(T.Rtopduration(ii) == durats);
            col = find(T.Rtoponset(ii) == onsets);  
            T.choiceArray(T.Rtrialtype(ii)+1) = T.choiceArray(T.Rtrialtype(ii)+1) + 1;
            T.choiceMatrix(row,col) = T.choiceMatrix(row,col) + 1;
        end
    %     T.choiceMatrix(row,col) = T.choiceMatrix(row,col) + 1;  % THIS NEEDS TO BE INSIDE VALID CHOICE WRAPPER

        T.combotype(ii) = find(ismember(combos,sort([T.Ltrialtype(ii),T.Rtrialtype(ii)]),'rows'));
        if ~isnan(T.choice(ii))
    %         T.combotype(ii) = find(ismember(combos,sort([T.Ltrialtype(ii),T.Rtrialtype(ii)]),'rows')); % THIS NEEDS TO BE OUTSIDE THE WRAPPER
            if T.Rtrialtype(ii) < T.Ltrialtype(ii)
                if T.choice(ii) == 1
                    T.choiceTypeMatrix(T.combotype(ii),2) = T.choiceTypeMatrix(T.combotype(ii),2) + 1;
                elseif T.choice(ii) == 2
                    T.choiceTypeMatrix(T.combotype(ii),1) = T.choiceTypeMatrix(T.combotype(ii),1) + 1;                    
                end
            else
                T.choiceTypeMatrix(T.combotype(ii),T.choice(ii)) = T.choiceTypeMatrix(T.combotype(ii),T.choice(ii)) + 1;
            end
        end
    end
    % T.p_choiceMatrix = T.choiceMatrix/T.trials; % THIS ISNT VALID, since T.trials assumes all valid trials
    T.p_choiceMatrix = T.choiceMatrix/sum(~isnan(T.choice));

    % Cant do this... because there arent 5 trial types for each...
    duratsOnly = sum(T.choiceMatrix,2);
    onsetsOnly = sum(T.choiceMatrix,1);


    %% Save results to Structure
    m{s} = T;

end