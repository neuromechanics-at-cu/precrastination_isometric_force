function create_trial_list(session,subj,followup)
%% Create two text files of trial combinations
% One for familiarization
% One for the choice session
%
% Session == 1: Hills
% Session == 2: Valleys


if session == 1
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
    % Duration of Effort (5 conditions) (s)
    durats = [1,3,5];

elseif session == 2
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
    % Duration of Effort (5 conditions) (s)
    durats = [1,3,5];
end


% Initialize big matrix
bigMat = [];

% Run through combinations
for jj = 1:length(durats)
    for ii = 1:size(onsets,2)
        if ~isnan(onsets(jj,ii))
                Xtemp = [0;
                         rampuptime2base;
                         onsets(jj,ii);
                         onsets(jj,ii) + rampuptime2top;
                         onsets(jj,ii) + rampuptime2top + durats(jj);
                         onsets(jj,ii) + rampuptime2top + durats(jj) + rampdowntime2base;
                         totalduration - rampdowntime2end;
                         totalduration];
                Ytemp = [0;
                         baselevel;
                         baselevel;
                         toplevel;
                         toplevel;
                         baselevel;
                         baselevel;
                         0];
%               % Plot Profiles
%             figure
%             plot(Xtemp,Ytemp)
%             xlim([0 totalduration])
%             ylim([0 1.1])
            bigMat = [bigMat;
                      totalduration,rampuptime2base,rampuptime2top,rampdowntime2base,rampdowntime2end,baselevel,toplevel,durats(jj),onsets(jj,ii)];
            
        end
    end
end

if ~exist('followup') || strcmp(followup,'no')
    
    % % Write to CSV - trailtypes
    % dlmwrite('trialtypes.txt',bigMat,'delimiter','\t')
    
    % Write to CSV - Trial order
    trialorder_array = [];
    % NumTimes Experienced
    ntimes = 3;
    % Randomize Order
    exposureOrder = randperm(size(bigMat,1)) - 1;
    % Ordered exposure
    for ii = 1:length(exposureOrder)
        trialorder_array = [trialorder_array; ones(ntimes,1)*exposureOrder(ii)];
    end
    % random selection
    % need to gaurantee 3 of each condition, 
    ntimes_rand = 1;
    randomtrials = randperm(ntimes_rand*size(bigMat,1));% gaurantee n occurences
    randomtrials = mod(randomtrials,size(bigMat,1));% make it range from 0 to size(bigMat,1) - 1
    
    trialorder_array = [trialorder_array; randomtrials'];
    
    % Add trialnumbers
    trialorder_array = [[1:length(trialorder_array)]',trialorder_array];
    % csvwrite('trialorder',trialorder_array)
    % dlmwrite('trialorders.txt',trialorder_array,'delimiter','\t')
    
    % Write to CSV - Choice comparison
    % guarantee each combos is seen 5 times
    
    % The code below looks at 5 of each R/L bias.. probably a bit much.  To
    % git rid of those, uncomment the line below:
    % combos = combvec([0:size(bigMat,1)-1],[0:size(bigMat,1)-1]);
    combos = nchoosek([0:size(bigMat,1)-1],2);
    
    list = [];
    for ii = 1:5
        list = [list,randperm(size(combos,1))];
    end
    randlistind = randperm(length(list));
    randomizedlist = list(randlistind);
    choicearray = combos(randomizedlist,:); 
    %Add realization that gaurantees each combo is realized once
    realize = zeros(size(randomizedlist));
    % realize(1:size(combos,1)) = ones(1,size(combos,1));
    % INSTEAD ONLY REALIZE HALF OF CHOICES (FOR FATIGUE)... roughly 10% of
    % trials for 14 choices...
    half = round(0.5*size(combos,1));
    realizehalf = ones(1,half);
    realizetempind = randperm(size(combos,1));
    realizetemp = zeros(1,size(combos,1));
    realizetemp(1:half) = realizehalf;
    realizetemp = realizetemp(realizetempind);
    
    
    realize(1:size(combos,1)) = realizetemp;
    realize = realize(randlistind);
    
    % % randomize right/left hand side
    % % this isnt "pseudorandom" in the sense that it doesnt gaurantee that each
    % % is seen an equal amount of times on each side.
    swap = round(rand(size(choicearray,1),1)); 
    for ii = 1:length(swap)
        if swap(ii) 
            choicearray(ii,:) = [choicearray(ii,2),choicearray(ii,1)];
        end
    end
    % 
    % dlmwrite('choiceorder.txt',[1:size(choicearray,1);choicearray']','delimiter','\t')
    
    
    %% EXPOSURE TRIAL LIST AND SPECS
    % Organization:
    % trialnumber/trialtype/timestypeseen/totalduration/rampuptime/rampdowntime/baselevel/toplevel/duration/onsettime
    completetriallist = [];
    count = zeros(size(bigMat,1),1);
    for ii = 1:size(trialorder_array)
        type = trialorder_array(ii,2)+1;
        count(type) = count(type)+1;% increment count
        completetriallist = [completetriallist;
                             trialorder_array(ii,:),count(type),bigMat(type,:)];                
    end
    % Save as txt file
    if session == 2
        dlmwrite(['exposuretriallist_',subj,'2.txt'],completetriallist,'delimiter','\t')
    elseif session == 1
        dlmwrite(['exposuretriallist_',subj,'.txt'],completetriallist,'delimiter','\t')
    end
    
    % CHOICE TASK AND SPECS
    % Organization:
    % trialnumber/left-type1/right-type2/timescomboseen/timestype1seen/timestype2seen/...
    % totalduration1/rampuptime1/rampdowntime1/baselevel1/toplevel1/duration1/onsettime1/...
    % totalduration2/rampuptime2/rampdowntime2/baselevel2/toplevel2/duration2/onsettime2
    completechoicelist = [];
    count_combo = zeros(size(combos,1),1);
    count_type1 = zeros(size(bigMat,1),1);
    count_type2 = zeros(size(bigMat,1),1);
    for ii = 1:size(choicearray)
        type1 = choicearray(ii,1)+1;
        type2 = choicearray(ii,2)+1;
        combotype = find(ismember(combos,sort(choicearray(ii,:)),'rows'));
        count_combo(combotype) = count_combo(combotype)+1;% increment count
        count_type1(type1) = count_type1(type1)+1;% increment count
        count_type2(type2) = count_type2(type2)+1;% increment count
        completechoicelist = [completechoicelist;
                             ii,choicearray(ii,:),count_combo(combotype),count_type1(type1),count_type2(type2),bigMat(type1,:),bigMat(type2,:),realize(ii)];                
    end
    % Save as txt file
    if session == 2
        dlmwrite(['choicetriallist_',subj,'2.txt'],completechoicelist,'delimiter','\t')
    elseif session == 1
        dlmwrite(['choicetriallist_',subj,'.txt'],completechoicelist,'delimiter','\t')        
    end

elseif strcmp(followup,'allcombosonce')
    combos = nchoosek([0:size(bigMat,1)-1],2);
    
%     list = [];
%     for ii = 1:5
%         list = [list,randperm(size(combos,1))];
%     end
    list = randperm(size(combos,1));
    randlistind = randperm(length(list));
    randomizedlist = list(randlistind);
    choicearray = combos(randomizedlist,:); 
    %Add realization that gaurantees each combo is realized once
    realize = zeros(size(randomizedlist));
    % realize(1:size(combos,1)) = ones(1,size(combos,1));
    % INSTEAD ONLY REALIZE HALF OF CHOICES (FOR FATIGUE)... roughly 10% of
    % trials for 14 choices...
%     half = round(0.5*size(combos,1));
    proportion = round(0.1*size(list,2));
    realizeproportion = ones(1,proportion);
    realizetempind = randperm(size(combos,1));
    realizetemp = zeros(1,size(combos,1));
    realizetemp(1:proportion) = realizeproportion;
    realizetemp = realizetemp(realizetempind);
   
    realize(1:size(combos,1)) = realizetemp;
    realize = realize(randlistind);
    
    % % randomize right/left hand side
    % % this isnt "pseudorandom" in the sense that it doesnt gaurantee that each
    % % is seen an equal amount of times on each side.
    swap = round(rand(size(choicearray,1),1)); 
    for ii = 1:length(swap)
        if swap(ii) 
            choicearray(ii,:) = [choicearray(ii,2),choicearray(ii,1)];
        end
    end

    completechoicelist = [];
    count_combo = zeros(size(combos,1),1);
    count_type1 = zeros(size(bigMat,1),1);
    count_type2 = zeros(size(bigMat,1),1);
    for ii = 1:size(choicearray)
        type1 = choicearray(ii,1)+1;
        type2 = choicearray(ii,2)+1;
        combotype = find(ismember(combos,sort(choicearray(ii,:)),'rows'));
        count_combo(combotype) = count_combo(combotype)+1;% increment count
        count_type1(type1) = count_type1(type1)+1;% increment count
        count_type2(type2) = count_type2(type2)+1;% increment count
        completechoicelist = [completechoicelist;
                             ii,choicearray(ii,:),count_combo(combotype),count_type1(type1),count_type2(type2),bigMat(type1,:),bigMat(type2,:),realize(ii)];                
    end    
    % Save as txt file
    if session == 2
        dlmwrite(['choicetriallist_',subj,'2_followup.txt'],completechoicelist,'delimiter','\t')
    elseif session == 1
        dlmwrite(['choicetriallist_',subj,'_followup.txt'],completechoicelist,'delimiter','\t')        
    end    
elseif strcmp(followup,'subset')
    combos = nchoosek([0:size(bigMat,1)-1],2);
    newchoices = [6,10,67,59,63,85,26,71,89,5,66,88,1,47]; % ORIGINAL
    newchoices = sort(newchoices);
    combos = combos(newchoices,:);
    
    list = [];
    for ii = 1:5
        list = [list,randperm(size(combos,1))];
    end
    randlistind = randperm(length(list));
    randomizedlist = list(randlistind);
    choicearray = combos(randomizedlist,:); 
    %Add realization that gaurantees each combo is realized once
    realize = zeros(size(randomizedlist));
    % realize(1:size(combos,1)) = ones(1,size(combos,1));
    % INSTEAD ONLY REALIZE HALF OF CHOICES (FOR FATIGUE)... roughly 10% of
    % trials for 14 choices...
    half = round(0.5*size(combos,1));
    realizehalf = ones(1,half);
    realizetempind = randperm(size(combos,1));
    realizetemp = zeros(1,size(combos,1));
    realizetemp(1:half) = realizehalf;
    realizetemp = realizetemp(realizetempind);
    
    
    realize(1:size(combos,1)) = realizetemp;
    realize = realize(randlistind);
    
    
    % % randomize right/left hand side
    % % this isnt "pseudorandom" in the sense that it doesnt gaurantee that each
    % % is seen an equal amount of times on each side.
    swap = round(rand(size(choicearray,1),1)); 
    for ii = 1:length(swap)
        if swap(ii) 
            choicearray(ii,:) = [choicearray(ii,2),choicearray(ii,1)];
        end
    end

    completechoicelist = [];
    count_combo = zeros(size(combos,1),1);
    count_type1 = zeros(size(bigMat,1),1);
    count_type2 = zeros(size(bigMat,1),1);
    for ii = 1:size(choicearray)
        type1 = choicearray(ii,1)+1;
        type2 = choicearray(ii,2)+1;
        combotype = find(ismember(combos,sort(choicearray(ii,:)),'rows'));
        count_combo(combotype) = count_combo(combotype)+1;% increment count
        count_type1(type1) = count_type1(type1)+1;% increment count
        count_type2(type2) = count_type2(type2)+1;% increment count
        completechoicelist = [completechoicelist;
                             ii,choicearray(ii,:),count_combo(combotype),count_type1(type1),count_type2(type2),bigMat(type1,:),bigMat(type2,:),realize(ii)];                
    end    
    % Save as txt file
    if session == 2
        dlmwrite(['choicetriallist_',subj,'2_followup.txt'],completechoicelist,'delimiter','\t')
    elseif session == 1
        dlmwrite(['choicetriallist_',subj,'_followup.txt'],completechoicelist,'delimiter','\t')        
    end    
end