%% Initialize Data and Parameters
% Exploring fits using full data set.


%% Load Data
load("../../Data/Processed Data/processed_data.mat")
addpath('../Supporting Functions')

%% Initialize Params
reward = 0;

% Utility Model Type Descriptions
% 1: No Temporal Sensitivity, Force Only
% 2: Temporal Sensitivity, Force Only
% 3: No temporal Sensitivity, Combined
% 4: Temporal Sensitivity, Combined
% 5: Temporal Sensitivity, Separate Discounting
% 6: Temporal Sensitivity of Task Only, Time Neutral Force
%%%NOTE: Comment this out if running from supercomputer.
utilitymodeltypes = 1;
%%%%%

% Task Model Type Descriptions
% 0 - Fdot
% 1 - Carry 
% 2 - Hill/Valley Only
taskmodeltypes = 0; 



% Discount Model Types
% 0 - Exponential
% 1 - Hyperbolic
% 2 - Log (TODO)
discountmodeltypes = 2; 

% Fit Search Parameters
% Cost/Minimization Type
% 1: RMS, 2: Likelihood
costtype = 2;
% options = optimset('TolFun',1e-6,'MaxIter',100000,'MaxFunEVals',100000);
options = optimset('TolFun',1e-6,'MaxIter',10000,'MaxFunEVals',10000);% update after stuck on subj6
nonlcon = [];    
AA = [];
bb = [];
Aeq = [];
beq = [];

% For Exponential
% gamma_lb = -3;
gamma_lb = -0.7; % 2023 Update based on SOL (and randomized starts causing NaN returns)
% gamma_ub = 3;
gamma_ub = 0.7;% 2023 New update for the peicewise discount fxn (2)
beta_lb = 0;
% beta_ub = 1; % Use this beta for exponential
beta_ub = 1e6; % Update 2023, normed

% Normed Approach
rel_weight_lb = 0;
rel_weight_ub = 1;

% For Hyperbolic
% gamma_lb = -1e3;
% gamma_ub = 1e3;
% beta_lb = 0;
% beta_ub = 1;
% rel_weight_lb = 1;
% rel_weight_ub = 1e3;

% Number of Initial Conditions 
init_gammaPhysIs = 1;
init_gammaTaskIs = 1;
init_betais = 1;
init_relweight = 1;

% Initialize Outputs
% sol = nan(length(utilitymodeltypes),length(taskmodeltypes),length(discountmodeltypes),...
%             length(costtype),length(subjs),init_gammaPhysIs,init_gammaTaskIs,...
%             init_betais,init_relweight,4);
% fval = nan(length(utilitymodeltypes),length(taskmodeltypes),length(discountmodeltypes),...
%             length(costtype),length(subjs),init_gammaPhysIs,init_gammaTaskIs,...
%             init_betais,init_relweight);
sol = nan(length(subjs),4);
fval = nan(length(subjs));

% Restarts
n_restarts = 10;

%% Iterate Through Model Fits
% for ii = 1:length(utilitymodeltypes)
%     for kk = 1:length(taskmodeltypes)
%         for qq = 1:length(discountmodeltypes)
%             for rr = 1:length(costtype)
ii = 1;kk=1;qq=1;rr=1;
    
% Initial Parameters for each Model Type
switch utilitymodeltypes(ii)
    case 1 % Time Neutral, Force Only
        ub = [0,0,beta_ub,0];
        lb = [0,0,beta_lb,0];
        gammaPhysIs=0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight=1;
    case 2 % Time Sensitive, Force Only
        ub = [gamma_ub,0,beta_ub,0];
        lb = [gamma_lb,0,beta_lb,0];
        gammaPhysIs=0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight=1;
    case 3 % Time Neutral, Force and Force Dot
        ub = [0,0,beta_ub,rel_weight_ub];
        lb = [0,0,beta_lb,rel_weight_lb];
        gammaPhysIs= 0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight = 1;
    case 4 % Time Sensitive Force and Force Dot are the same
        ub = [gamma_ub,0,beta_ub,rel_weight_ub];
        lb = [gamma_lb,0,beta_lb,rel_weight_lb];
        gammaPhysIs= 0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight = 1;
    case 5 % Time Sensitive, Force and Force Dot are independent
        ub = [gamma_ub,gamma_ub,beta_ub,rel_weight_ub];
        lb = [gamma_lb,gamma_lb,beta_lb,rel_weight_lb];
        gammaPhysIs= 0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight = 1;
    case 6 % Time Neutral Force, Time Sensitive Force Dot
        ub = [0,gamma_ub,beta_ub,rel_weight_ub];
        lb = [0,gamma_lb,beta_lb,rel_weight_lb];
        gammaPhysIs= 0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight = 1;
    case 7 % Minimize Time in State
        ub = [0,0,beta_ub,-1];
        lb = [0,gamma_lb,beta_lb,-1];
        gammaPhysIs= 0; % OVERALL TIME
        gammaTaskIs= -1e-3; %TIME IN TASK
        betais=1e-3;
        rel_weight = -1;
    case 8 % Minimize Time in State + Time Sensitivity 
        ub = [gamma_ub,0,beta_ub,rel_weight_ub];
        lb = [gamma_lb,gamma_lb,beta_lb,rel_weight_lb];
        gammaPhysIs= 0; % OVERALL TIME
        gammaTaskIs=-1e-3; %TIME IN TASK
        betais=1e-3;
        rel_weight = 1;
    case 9 % Time Sensitive Force Dot Only
        ub = [0,gamma_ub,beta_ub,1];
        lb = [0,gamma_lb,beta_lb,1];
        gammaPhysIs= 0;
        gammaTaskIs=0;
        betais=1e-3;
        rel_weight = 1;
end

parfor jj = 1:length(subjs)    
%     for mm = 1:length(gammaPhysIs)
%         for pp = 1:length(gammaTaskIs)
%             for nn = 1:length(betais)
%                 for ll = 1:length(rel_weight)
    mm=1;pp=1;nn=1;ll=1;

    best_fval = 1e4;
    disp(['Subj ',num2str(jj),', Model ',num2str(utilitymodeltypes(ii)),...
    ', Task Model ',num2str(taskmodeltypes(kk)),', Discount Type ',num2str(discountmodeltypes(qq)),...
    ', Cost Type ', num2str(costtype(rr)),...
    ', Gamma Physical #',num2str(mm),' of ',num2str(length(gammaPhysIs)),...
    ', Gamma Task #',num2str(pp),' of ',num2str(length(gammaTaskIs)),...
    ' and Beta #',num2str(nn),' of ',num2str(length(nn)),...
      ' and relweight #',num2str(ll),' of ',num2str(length(ll))])
    for zz = 1:n_restarts
        disp(['Random Start #',num2str(zz),'Subj ',num2str(jj)])

         x0 = [(ub(1)-lb(1)).*rand(1,1) + lb(1),...
            (ub(2)-lb(2)).*rand(1,1) + lb(2),...
            1e-3,...
            (ub(4)-lb(4)).*rand(1,1) + lb(4)]';

        while isnan(fitSubjectiveUtilityModel_normed(x0,Mchoice{jj},reward,utilitymodeltypes(ii),costtype(rr),taskmodeltypes(kk),discountmodeltypes(qq)))
            x0 = [(ub(1)-lb(1)).*rand(1,1) + lb(1),...
                (ub(2)-lb(2)).*rand(1,1) + lb(2),...
                1e-3,...
                (ub(4)-lb(4)).*rand(1,1) + lb(4)]';
            disp('Choosing new start point')
        end

        if length(Mchoice) == 1
            [sol_temp,fval_temp,exitflag] = ...
                fmincon(@(x) fitSubjectiveUtilityModel_normed(x,Mchoice,reward,utilitymodeltypes(ii),costtype(rr),taskmodeltypes(kk),discountmodeltypes(qq)),...
                x0,AA,bb,Aeq,beq,lb,ub,nonlcon,options);                                        
        else
            [sol_temp,fval_temp,exitflag] = ...
                fmincon(@(x) fitSubjectiveUtilityModel_normed(x,Mchoice{jj},reward,utilitymodeltypes(ii),costtype(rr),taskmodeltypes(kk),discountmodeltypes(qq)),...
                x0,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
        end
        if exitflag < 1
            disp('Bad Exit.')
        end
        % Keep Best Solution
        if fval_temp < best_fval
%             sol(ii,kk,qq,rr,jj,mm,pp,nn,ll,:) = sol_temp;
            sol(jj,:) = sol_temp;
%             fval(ii,kk,qq,rr,jj,mm,pp,nn,ll) = fval_temp;
            fval(jj) = fval_temp;
            best_fval = fval_temp;
            disp(['New Best Found,',num2str(zz)])
        end
    end
end

save(['model_',num2str(utilitymodeltypes),'.mat'])
