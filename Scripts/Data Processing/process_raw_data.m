%% Process Raw Data
% Choose Subset 
% 0: all
% 1: last two choices
% 2: first two choices
subset = 0;

%% List of Subjects and Test Session Conditions
% Ages, dates, handedness omitted

% List only subjects who have completed both sessions
subjs = {'WLH','KRT','DNF','MMM','MAH','MEE','ALS','WAS','AFT','DAM','BLM','CAL','JEL','DUF','LOG',...
    'CRG','PKJ','ARR','JFF','RAP','CNK','SHG','JKB','VNN','GKH'};
% Session "A", Hills, MVC
mvc1A = [2.485,2.158,2.803,1.913,1.373,1.814,1.034,3.361,2.481,4.440,2.363,2.601,3.490,2.324,2.492,...
    4.990,2.925,2.004,2.802,2.395,1.814,2.376,1.095,2.412,1.845];
mvc2A = [2.550,2.204,2.847,2.021,1.616,2.384,0.851,3.728,2.3333,4.602,2.041,2.709,3.290,2.541,2.597,...
    3.875,3.351,1.927,2.817,2.480,1.943,2.360,1.237,2.221,2.229];
mvc3A = [2.580,2.239,3.089,2.217,1.751,2.612,1.000,3.651,2.801,4.969,2.629,3.418,2.424,2.927,2.671,...
    3.804,3.095,2.003,3.173,2.591,2.076,2.393,1.310,2.409,2.276];
% Session "A", Hills, Subjective Fatigue
fatigue1A = [1,2,1.5,1,1,2,1,1,1,2.5,1,2,1,2,1,...
    1,1,1,2,1,1,1,1,1,1];
fatigue2A = [3,4,3,3,2,3,3,2.5,2,3,2,3,3,4,3,...
    2,2,1.5,2,2,4.5,1,3,3,2.5];
fatigue3A = [3.5,3,3.5,2,2,3,3.5,2.5,3,3,2,3,4.5,3,2,...
    2.5,3,2,3,2,3,1,3.5,4,3];
% Session "A" Elbow Angle
el_angA = [233.7,234.3,234.67,234.481,232.690,232.993,235.959,232.971,235.827,235.997,235.871,236.025,235.750,235.415,236.124,...
    233.080,232.993,233.322,233.366,233.509,233.168,232.498,233.179,233.251,233.097];

% Session "B", Valleys, MVC
mvc1B = [2.373,2.247,2.878,2.332,2.364,2.875,1.251,3.591,2.542,2.875,2.157,2.749,3.614,2.670,2.174,...
    4.293,3.283,1.794,2.661,2.568,2.035,1.942,1.593,1.988,2.327];
mvc2B = [2.443,2.209,2.739,2.431,2.451,3.221,1.236,3.329,2.286,3.221,2.082,2.546,1.969,2.938,2.220,...
    4.762,3.247,2.000,2.779,2.372,1.993,2.045,1.398,1.994,2.442];
mvc3B = [2.475,2.209,2.860,2.472,1.751,2.343,0.977,3.605,2.801,3.861,2.368,2.711,NaN,2.965,2.181,...
    4.611,3.683,1.881,2.963,2.312,2.172,1.984,1.761,2.053,2.772];
% Session "B", Valleys, Subjective Fatigue
fatigue1B = [1,1,1,1,1,1.5,1,1,1,2,1,2,1,2,2,...
    1,1,1,1,1,1,1,1,1,1];
fatigue2B = [2,3,2,2,2,2,3,2,2,2,2,3,4.5,3,3,...
    2.5,3,2,3,3,3,2,2,3,2];
fatigue3B = [3,3,3,2,2,2,2,2.5,3,3,2,4,NaN,3,2,...
    3.5,3,1.5,3,2,2.5,1,2.5,4,2.5];
% Session "B" Elbow Angle
el_angB = [234.586,234.624,234.904,232.696,232.927,232.910,233.476,235.992,233.569,232.910,235.981,235.948,235.920,235.992,235.503,...
    233.443,233.767,233.756,233.086,233.168,233.306,233.361,233.256,233.179,233.015];

% Add Path to Supporting Functions
addpath('../../Scripts/Supporting Functions')

% Load Familiarization trials
M_A = processData_Familiarization('../../Data/Raw Data/Familiarization/Hill',subjs,mvc1A,el_angA);
M_B = processData_Familiarization('../../Data/Raw Data/Familiarization/Valley',subjs,mvc1B,el_angB);

% Load Choice trials
Mchoice_A = processData_Choice('../../Data/Raw Data/Choice/Hill',subjs,mvc2A,el_angA);
Mchoice_B = processData_Choice('../../Data/Raw Data/Choice/Valley',subjs,mvc2B,el_angB);

% Combine Data Structures
% [M,Mchoice] = combineDataStructs(M_A,M_B,Mchoice_A,Mchoice_B);
[M,Mchoice] = combineDataStructs_subset(M_A,M_B,Mchoice_A,Mchoice_B,subset);


% Save Data Struct
clear M_A M_B Mchoice_A Mchoice_B
save('../../Data/Processed Data/processed_from_raw.mat','-v7.3')