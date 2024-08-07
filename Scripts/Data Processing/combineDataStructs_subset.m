function [M,Mchoice] = combineDataStructs_subset(M_A,M_B,Mchoice_A,Mchoice_B,subsettype)

if subsettype == 1
    % Last two choices in session
    startind = 274;
    endind = 455;
    start_fam = 1;
    end_fam = 56;
elseif subsettype == 2
    % First two choices in session
    startind = 1;
    endind = 182;
    start_fam = 1;
    end_fam = 56;    
else
    startind = 1;
    endind = 455;
    start_fam = 1;
    end_fam = 56;
end

if size(M_A,2) ~= size(M_B,2) || size(Mchoice_A,2) ~= size(Mchoice_B,2)
    error('Data structures are not the same size.')
end


for ii = 1:size(M_A,2)
    
    % Familiarization 
%     M{ii}.trials = M_A{ii}.trials + M_B{ii}.trials;
    M{ii}.trials = 2*(endind - startind+1);

    M{ii}.framedata = [M_A{ii}.framedata(start_fam:end_fam)';
                       M_B{ii}.framedata(start_fam:end_fam)']';
    M{ii}.framedata = M{ii}.framedata(start_fam:end_fam)';
    M{ii}.trialnumber = [1:length(start_fam:end_fam),(length(start_fam:end_fam) + 1):(2*length(start_fam:end_fam))]';
    M{ii}.trialnumber_bysession = [1:length(start_fam:end_fam),1:length(start_fam:end_fam)]';
    M{ii}.trialtype = [M_A{ii}.trialtype(start_fam:end_fam);
                       max(M_A{ii}.trialtype(start_fam:end_fam)) + 1 + M_B{ii}.trialtype(start_fam:end_fam)]; 
    M{ii}.trialtypecount = [M_A{ii}.trialtypecount(start_fam:end_fam);
                            M_B{ii}.trialtypecount(start_fam:end_fam)];
    M{ii}.trialduration = [M_A{ii}.trialduration(start_fam:end_fam);
                           M_B{ii}.trialduration(start_fam:end_fam)];
    M{ii}.rampuptime2base = [M_A{ii}.rampuptime2base(start_fam:end_fam);
                             M_B{ii}.rampuptime2base(start_fam:end_fam)];
    M{ii}.rampuptime2top = [M_A{ii}.rampuptime2top(start_fam:end_fam);
                            M_B{ii}.rampuptime2top(start_fam:end_fam)];
    M{ii}.rampdowntime2base = [M_A{ii}.rampdowntime2base(start_fam:end_fam);
                               M_B{ii}.rampdowntime2base(start_fam:end_fam)];                      
    M{ii}.rampdowntime2end = [M_A{ii}.rampdowntime2end(start_fam:end_fam);
                               M_B{ii}.rampdowntime2end(start_fam:end_fam)];     
    M{ii}.baselevel = [M_A{ii}.baselevel(start_fam:end_fam);
                       M_B{ii}.baselevel(start_fam:end_fam)];   
    M{ii}.toplevel = [M_A{ii}.toplevel(start_fam:end_fam);
                       M_B{ii}.toplevel(start_fam:end_fam)];    
    M{ii}.topduration = [M_A{ii}.topduration(start_fam:end_fam);
                       M_B{ii}.topduration(start_fam:end_fam)];
    M{ii}.toponset = [M_A{ii}.toponset(start_fam:end_fam);
                       M_B{ii}.toponset(start_fam:end_fam)];
    M{ii}.accuracy = [M_A{ii}.accuracy(start_fam:end_fam);
                      M_B{ii}.accuracy(start_fam:end_fam)];   
    M{ii}.peakRateOfForceUphill_px = [M_A{ii}.peakRateOfForceUphill_px(start_fam:end_fam);
                                      M_B{ii}.peakRateOfForceUphill_px(start_fam:end_fam)];
    M{ii}.peakRateOfForceDownhill_px = [M_A{ii}.peakRateOfForceDownhill_px(start_fam:end_fam);
                                      M_B{ii}.peakRateOfForceDownhill_px(start_fam:end_fam)];
    M{ii}.peakRateOfForceUphill_N = [M_A{ii}.peakRateOfForceUphill_N(start_fam:end_fam);
                                      M_B{ii}.peakRateOfForceUphill_N(start_fam:end_fam)];
    M{ii}.peakRateOfForceDownhill_N = [M_A{ii}.peakRateOfForceDownhill_N(start_fam:end_fam);
                                      M_B{ii}.peakRateOfForceDownhill_N(start_fam:end_fam)];
    M{ii}.failedTrials = [M_A{ii}.failedTrials;
                          M_A{ii}.trials+M_B{ii}.failedTrials];
    % Accuracy Metrics
    M{ii}.mean_accuracytot = [M_A{ii}.mean_accuracytot,M_B{ii}.mean_accuracytot]';
    M{ii}.mean_accuracytot_last1 = [M_A{ii}.mean_accuracytot_last1,M_B{ii}.mean_accuracytot_last1]';
    M{ii}.mean_accuracytot_lastrand = [M_A{ii}.mean_accuracytot_lastrand,M_B{ii}.mean_accuracytot_lastrand]';
    M{ii}.forcetimeintegral_trialtype = [M_A{ii}.forcetimeintegral_trialtype,M_B{ii}.forcetimeintegral_trialtype]';

    % Choice
    Mchoice{ii}.trials = 2*(endind - startind+1);
    Mchoice{ii}.framedata = [Mchoice_A{ii}.framedata(startind:endind)';
                             Mchoice_B{ii}.framedata(startind:endind)']';
    Mchoice{ii}.framedata = Mchoice{ii}.framedata';
    Mchoice{ii}.trialnumber = [1:length(startind:endind),(length(startind:endind) + 1):(2*length(startind:endind))]';
    Mchoice{ii}.trialnumber_bysession = [1:length(startind:endind),1:length(startind:endind)]';
    Mchoice{ii}.Ltrialtype = [Mchoice_A{ii}.Ltrialtype(startind:endind);
                              max(Mchoice_A{ii}.Ltrialtype(startind:endind)) + 1 + Mchoice_B{ii}.Ltrialtype(startind:endind)];
    Mchoice{ii}.Rtrialtype = [Mchoice_A{ii}.Rtrialtype(startind:endind);
                              max(Mchoice_A{ii}.Rtrialtype(startind:endind)) + 1 + Mchoice_B{ii}.Rtrialtype(startind:endind)];
    Mchoice{ii}.choicetypecount = [Mchoice_A{ii}.choicetypecount(startind:endind);
                                   Mchoice_B{ii}.choicetypecount(startind:endind)];
    Mchoice{ii}.Ltrialtypecount = [Mchoice_A{ii}.Ltrialtypecount(startind:endind);
                                   Mchoice_B{ii}.Ltrialtypecount(startind:endind)];                               
    Mchoice{ii}.Rtrialtypecount = [Mchoice_A{ii}.Rtrialtypecount(startind:endind);
                                   Mchoice_B{ii}.Rtrialtypecount(startind:endind)]; 
    Mchoice{ii}.Ltrialduration = [Mchoice_A{ii}.Ltrialduration(startind:endind);
                           Mchoice_B{ii}.Ltrialduration(startind:endind)];
    Mchoice{ii}.Lrampuptime2base = [Mchoice_A{ii}.Lrampuptime2base(startind:endind);
                             Mchoice_B{ii}.Lrampuptime2base(startind:endind)];
    Mchoice{ii}.Lrampuptime2top = [Mchoice_A{ii}.Lrampuptime2top(startind:endind);
                            Mchoice_B{ii}.Lrampuptime2top(startind:endind)];
    Mchoice{ii}.Lrampdowntime2base = [Mchoice_A{ii}.Lrampdowntime2base(startind:endind);
                               Mchoice_B{ii}.Lrampdowntime2base(startind:endind)];                      
    Mchoice{ii}.Lrampdowntime2end = [Mchoice_A{ii}.Lrampdowntime2end(startind:endind);
                               Mchoice_B{ii}.Lrampdowntime2end(startind:endind)];     
    Mchoice{ii}.Lbaselevel = [Mchoice_A{ii}.Lbaselevel(startind:endind);
                            Mchoice_B{ii}.Lbaselevel(startind:endind)];   
    Mchoice{ii}.Ltoplevel = [Mchoice_A{ii}.Ltoplevel(startind:endind);
                            Mchoice_B{ii}.Ltoplevel(startind:endind)];    
    Mchoice{ii}.Ltopduration = [Mchoice_A{ii}.Ltopduration(startind:endind);
                            Mchoice_B{ii}.Ltopduration(startind:endind)];
    Mchoice{ii}.Ltoponset = [Mchoice_A{ii}.Ltoponset(startind:endind);
                            Mchoice_B{ii}.Ltoponset(startind:endind)];
    Mchoice{ii}.Rtrialduration = [Mchoice_A{ii}.Rtrialduration(startind:endind);
                                Mchoice_B{ii}.Rtrialduration(startind:endind)];
    Mchoice{ii}.Rrampuptime2base = [Mchoice_A{ii}.Rrampuptime2base(startind:endind);
                                    Mchoice_B{ii}.Rrampuptime2base(startind:endind)];
    Mchoice{ii}.Rrampuptime2top = [Mchoice_A{ii}.Rrampuptime2top(startind:endind);
                                    Mchoice_B{ii}.Rrampuptime2top(startind:endind)];
    Mchoice{ii}.Rrampdowntime2base = [Mchoice_A{ii}.Rrampdowntime2base(startind:endind);
                                    Mchoice_B{ii}.Rrampdowntime2base(startind:endind)];                      
    Mchoice{ii}.Rrampdowntime2end = [Mchoice_A{ii}.Rrampdowntime2end(startind:endind);
                                    Mchoice_B{ii}.Rrampdowntime2end(startind:endind)];     
    Mchoice{ii}.Rbaselevel = [Mchoice_A{ii}.Rbaselevel(startind:endind);
                            Mchoice_B{ii}.Rbaselevel(startind:endind)];   
    Mchoice{ii}.Rtoplevel = [Mchoice_A{ii}.Rtoplevel(startind:endind);
                            Mchoice_B{ii}.Rtoplevel(startind:endind)];    
    Mchoice{ii}.Rtopduration = [Mchoice_A{ii}.Rtopduration(startind:endind);
                            Mchoice_B{ii}.Rtopduration(startind:endind)];
    Mchoice{ii}.Rtoponset = [Mchoice_A{ii}.Rtoponset(startind:endind);
                            Mchoice_B{ii}.Rtoponset(startind:endind)];
    Mchoice{ii}.realization = [Mchoice_A{ii}.realization(startind:endind);
                               Mchoice_B{ii}.realization(startind:endind)];
    Mchoice{ii}.choice = [Mchoice_A{ii}.choice(startind:endind);
                          Mchoice_B{ii}.choice(startind:endind)];
    Mchoice{ii}.choicetime = [Mchoice_A{ii}.choicetime(startind:endind);
                              Mchoice_B{ii}.choicetime(startind:endind)];   
    Mchoice{ii}.accuracy = [Mchoice_A{ii}.accuracy(startind:endind);
                            Mchoice_B{ii}.accuracy(startind:endind)];
    Mchoice{ii}.combotype = [Mchoice_A{ii}.combotype(startind:endind)';
                             max(Mchoice_A{ii}.combotype(startind:endind)) + Mchoice_B{ii}.combotype(startind:endind)']';
    Mchoice{ii}.combotype = Mchoice{ii}.combotype';
                         
    Mchoice{ii}.choiceTypeMatrix = [Mchoice_A{ii}.choiceTypeMatrix;
                                    Mchoice_B{ii}.choiceTypeMatrix];
    Mchoice{ii}.failedTrials = [Mchoice_A{ii}.failedTrials;
                          Mchoice_A{ii}.trials+Mchoice_B{ii}.failedTrials];
    Mchoice{ii}.Lforcetimeintegral = [Mchoice_A{ii}.Lforcetimeintegral(startind:endind),...
                                        Mchoice_B{ii}.Lforcetimeintegral(startind:endind)]';
    Mchoice{ii}.Rforcetimeintegral = [Mchoice_A{ii}.Rforcetimeintegral(startind:endind),...
                                        Mchoice_B{ii}.Rforcetimeintegral(startind:endind)]';
    Mchoice{ii}.max_dxpx_dt = [Mchoice_A{ii}.max_dxpx_dt(startind:endind),...
                                        Mchoice_B{ii}.max_dxpx_dt(startind:endind)]';
    Mchoice{ii}.delibtime = [Mchoice_A{ii}.delibtime(startind:endind),...
                             Mchoice_B{ii}.delibtime(startind:endind)]';
    Mchoice{ii}.delibtime_normed = [Mchoice_A{ii}.delibtime_normed(startind:endind),...
                             Mchoice_B{ii}.delibtime_normed(startind:endind)]';
    Mchoice{ii}.delibtime_nooutliers = [Mchoice_A{ii}.delibtime_nooutliers(startind:endind),...
                             Mchoice_B{ii}.delibtime_nooutliers(startind:endind)]';
    % This is normed by all choices, not the subset specified above
    Mchoice{ii}.delibtime_nooutliers_normed = [Mchoice_A{ii}.delibtime_nooutliers_normed(startind:endind),...
                             Mchoice_B{ii}.delibtime_nooutliers_normed(startind:endind)]';
    % New From Revision - Familiarization Trial Performance
    Mchoice{ii}.avg_force_base1 = [Mchoice_A{ii}.avg_force_base1(startind:endind),Mchoice_B{ii}.avg_force_base1(startind:endind)]';
    Mchoice{ii}.std_force_base1 = [Mchoice_A{ii}.std_force_base1(startind:endind),Mchoice_B{ii}.std_force_base1(startind:endind)]';
    Mchoice{ii}.avg_force_hillvalley = [Mchoice_A{ii}.avg_force_hillvalley(startind:endind),Mchoice_B{ii}.avg_force_hillvalley(startind:endind)]';
    Mchoice{ii}.std_force_hillvalley = [Mchoice_A{ii}.std_force_hillvalley(startind:endind),Mchoice_B{ii}.std_force_hillvalley(startind:endind)]';
    Mchoice{ii}.avg_force_base2 = [Mchoice_A{ii}.avg_force_base2(startind:endind),Mchoice_B{ii}.avg_force_base2(startind:endind)]';
    Mchoice{ii}.std_force_base2 = [Mchoice_A{ii}.std_force_base2(startind:endind),Mchoice_B{ii}.std_force_base2(startind:endind)]';
    Mchoice{ii}.overshoot_ramp1 = [Mchoice_A{ii}.overshoot_ramp1(startind:endind),Mchoice_B{ii}.overshoot_ramp1(startind:endind)]';
    Mchoice{ii}.peakrate_mvc_ramp1 = [Mchoice_A{ii}.peakrate_mvc_ramp1(startind:endind),Mchoice_B{ii}.peakrate_mvc_ramp1(startind:endind)]';
    Mchoice{ii}.overshoot_ramp2 = [Mchoice_A{ii}.overshoot_ramp2(startind:endind),Mchoice_B{ii}.overshoot_ramp2(startind:endind)]';
    Mchoice{ii}.peakrate_mvc_ramp2 = [Mchoice_A{ii}.peakrate_mvc_ramp2(startind:endind),Mchoice_B{ii}.peakrate_mvc_ramp2(startind:endind)]';
end