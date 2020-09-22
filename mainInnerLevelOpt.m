%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code solves the muscle-force sharing problem using only the inner
% level optimization. Therefore the optimized muscle-tendon parameters need
% to be provided. Parameters obtained in one of the 9 two-level 
% optimization Problems can be chosen.
% The details of the study can be found in the following article:
% Serrancolí, G; Kinney, A; Fregly, B.J. Influence of Musculoskeletal Model 
% Parameter Values on Prediction of Accurate Knee Contact Forces during 
% Walking. Medial Engineering & Physics 2020. In Press. 
% https://doi.org/10.1016/j.medengphy.2020.09.004

% Author: Gil Serrancolí

% The user can change several options at the Options structure. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Folders
foldermainroot=pwd;
folderExpData=[foldermainroot '\ProcessedExperimentalData\'];

%% Specify the Options of the algorithm
Options.parentCalibration='P7';                     %Parent Calibration to load the outer level variables
Options.patterncycle='ngait_og4';                   %Name of the cycle to estimate muscle forces
Options.MuscleParam='Arnold';                      %Whether use muscle parameters values from Delp 1990 ('Delp') or Arnold 2010 ('Arnold')
Options.weightingminA_FM0=0;                        %Whether use FM0_maxi/FM0 as weighting factors to minimize activations;
Options.InnerActivationsp05=0;                      %Whether minimizing (a_model+0.5)^2 in the Inner Level Cost Function instead of a_model^2
Options.UBactivations=1;                            %Whether use muscle activation bounds into Solve for Muscle Activations algorithm
Options.match1KneeIDload=0;                         %0--> Match 3 Knee ID loads (flexion, adduction, inf-sup); 1--> Match 1 ID load (flexion)
Options.StartFromlM0Reas=0;                         %Whether start from approx. good initial conditions for lM0 and lTs
Options.MaxlMtilda1_6trials=0;                      %Whether use constant lM0 and lTs values. Values optimized fixing max(lMtilda)=1;
Options.MaxlMtilda1_3trials=0;                      %Whether use constant lM0 and lTs values. Values optimized fixing max(lMtilda)=1, only using ngait trial 1, 2 and 3;
Options.muscleParamScaled=1;                        %Whether use scaled lM0 and lTs values. Values obtained after scaling gait2392 (with Arnold lM0 and lTs).
Options.folderExpData=folderExpData;                %Folder where the experimental data are located

%% Muscle names

muscles={'addbrev' 'addlong' 'addmagProx' 'addmagMid' 'addmagDist' 'addmagIsch' 'bflh' 'bfsh' 'edl' 'ehl' 'fdl' 'fhl' 'gaslat' 'gasmed' 'gem' 'glmax1' 'glmax2' 'glmax3' 'glmed1' 'glmed2' 'glmed3' 'glmin1' 'glmin2' 'glmin3' 'grac' 'iliacus' 'pect' 'perbrev' 'perlong' 'pertert' 'piri' 'psoas' 'quadfem' 'rf' 'sart' 'semimem' 'semiten' 'soleus' 'tfl' 'tibant' 'tibpost' 'vasint' 'vaslat' 'vasmed'};
nmuscles=length(muscles);

%% Load literature muscle parameters
if strcmp(Options.MuscleParam,'Delp')
    lM0lit=[0.1330 0.1380 0.1450 0.1310 0.1210 0.0870 0.1090 0.1730 0.1020 0.1110 0.0340 0.0430 0.0640 0.0450 0.0240 0.1420 0.1470 0.1440 0.0535 0.0845 0.0646 0.0680 0.0560 0.0380 0.3520 0.1000 0.1330 0.0500 0.0490 0.0790 0.0260 0.1040 0.0540 0.0840 0.5790 0.0800 0.2010 0.0300 0.0950 0.0980 0.0310 0.0870 0.0840 0.0890];
    lTslit=[0.0200 0.1100 0.1000 0.2600 0.1300 0.0600 0.3410 0.1000 0.3450 0.3050 0.4000 0.3800 0.3850 0.4080 0.0390 0.1250 0.1270 0.1450 0.0780 0.0530 0.0530 0.0160 0.0260 0.0510 0.1400 0.0900 0.0010 0.1610 0.3450 0.1000 0.1150 0.1300 0.0240 0.3460 0.0400 0.3590 0.2620 0.2680 0.4250 0.2230 0.3100 0.1360 0.1570 0.1260]; 
    FM0lit=[286 418 578 444 312 346 717 402 341 108 310 322 448 1113 109 382 546 368 546 382 435 180 190 215 108 429 177 348 754 90 296 371 254 779 104 1030 328 2839 155 603 1270 1365 1871 1294];
    PennationAngles=[0 6 5 5 3 5 0 23 8 6 7 10 8 17 0 5 0 5 8 0 19 10 0 1 3 7 0 5 10 13 10 8 0 5 0 15 5 25 3 5 12 3 5 5]*pi/180;
elseif strcmp(Options.MuscleParam,'Arnold')
    lM0lit=[0.1031 0.1082 0.1056 0.1377 0.1772 0.1562 0.0976 0.1103 0.0693 0.0748 0.0446 0.0527 0.0588 0.051 0.024 0.1473 0.1569 0.1665 0.0733 0.0733 0.0733 0.068 0.056 0.038 0.2278 0.1066 0.133 0.0454 0.0508 0.079 0.026 0.1169 0.054 0.0759 0.403 0.069 0.193 0.044 0.095 0.0683 0.0378 0.0993 0.0994 0.0968];
    lTslit=[0.036 0.13 0.043 0.048 0.09 0.221 0.322 0.104 0.3673 0.3315 0.3777 0.356 0.382 0.4008 0.039 0.05 0.0733 0.0702 0.0565 0.066 0.046 0.016 0.026 0.051 0.169 0.094 0.001 0.1481 0.333 0.1 0.115 0.097 0.024 0.346 0.11 0.378 0.245 0.2815 0.45 0.241 0.2818 0.106 0.13 0.112];
    FM0lit=[303.7 399.5 324.2 324.2 324.2 324.2 705.2 315.8 345.4 165 274.4 436.8 606.4 1308 109 546.1 780.5 526.1 881.1 616.5 702 180 190 215 137.3 621.9 177 305.9 653.3 90 296 479.7 254 848.8 113.5 1162.7 301.9 3585.9 155 673.7 905.6 1024.2 2255.4 1443.7];
    PennationAngles=[0.10646508 0.12356931 0.38728856 0.25673793 0.2412045 0.20734512 0.20210913 0.2151991 0.18901916 0.16475908 0.23806291 0.29478611 0.21013764 0.17243853 0.0 0.38292524 0.38292524 0.38292524 0.3572689 0.3572689 0.3572689 0.17453293 0.0 0.01745329 0.14241887 0.24940755 0.0 0.20001473 0.24574236 0.2268928 0.17453293 0.1860521 0.0 0.24312436 0.02321288 0.26337018 0.22444934 0.49305551 0.05235988 0.16685348 0.23928464 0.07923795 0.32079152 0.51679199];
end
lTslit(27)=0.01;

%% Load and Calculate data

% Calculate Normalized Parameter Values for use with the Rigid Tendon Model
params = zeros(34,nmuscles);

    [b,c,d,ev,eF,f,g,h] = GetNormalizedParameterValues('exp');
    for m = 1:nmuscles
        params(6:40,m) = [b; c; d; ev; eF; f; g; h];
    end

params(1,:) = FM0lit;params(1,:)=2*params(1,:);
params(2,:) = lM0lit;
params(3,:) = lTslit;
params(4,:) = PennationAngles*180/pi; % in degrees
params(5,:) = 10*lM0lit; %vMax in m/s (vMax0=10);
FM0=2*params(1,:);

switch Options.parentCalibration
    case 'P1'
        S=load('Results_OuterLevel/Solution_P1_varyinglM0.mat');
    case 'P2'
        S=load('Results_OuterLevel/Solution_P2_varyinglTs.mat');
    case 'P3'
        S=load('Results_OuterLevel/Solution_P3_varyingma.mat');
    case 'P4'
        S=load('Results_OuterLevel/Solution_P4_varyinglM0lTs.mat');
    case 'P5'
        S=load('Results_OuterLevel/Solution_P5_varyinglM0ma.mat');
    case 'P6'
        S=load('Results_OuterLevel/Solution_P6_varyinglTsma.mat');
    case 'P7'
        S=load('Results_OuterLevel/Solution_P7_varyinglM0lTsma.mat');
end
params_mod=S.params_mod;
patterncycle=Options.patterncycle;

% Load Inverse Dynamics Loads
cd(folderExpData);
if strcmp(patterncycle,'ngait_og4')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og4.mat');
elseif strcmp(patterncycle,'ngait_og5')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og5.mat');
elseif strcmp(patterncycle,'ngait_og7')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og7.mat');
elseif strcmp(patterncycle,'ngait_og1')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og1.mat');
elseif strcmp(patterncycle,'ngait_og2')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og2.mat');
elseif strcmp(patterncycle,'ngait_og3')
    ID_noCF=load('input_ID_withoutContactForces_ngait_og3.mat');
elseif strcmp(patterncycle,'2legsquat1')
    ID_noCF=load('input_ID_withoutContactForces_2legsquat1.mat');
elseif strcmp(patterncycle,'chairrise1')
    ID_noCF=load('input_ID_withoutContactForces_chairrise1.mat');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial1')
    ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial1.mat');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial2')
    ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial2.mat');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial3')
    ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial3.mat');
elseif strcmp(patterncycle,'openfe1')
    ID_noCF=load('input_ID_withoutContactForces_openfe1.mat');
elseif strcmp(patterncycle,'bouncy1')
    ID_noCF=load('input_ID_withoutContactForces_bouncy1.mat');
elseif strcmp(patterncycle,'medthrust2')
    ID_noCF=load('input_ID_withoutContactForces_medthrust2.mat');
elseif strcmp(patterncycle,'moderatecrouch3')
    ID_noCF=load('input_ID_withoutContactForces_moderatecrouch3.mat');
else
    keyboard;
end

%% Load IK data
if strcmp(patterncycle,'ngait_og4')
    IKdata=importdata('IK_ngait_og4_with_fluoro.mot');
elseif strcmp(patterncycle,'ngait_og5')
    IKdata=importdata('IK_ngait_og5_with_fluoro.mot');
elseif strcmp(patterncycle,'ngait_og7')
    IKdata=importdata('IK_ngait_og7_with_fluoro.mot');
elseif strcmp(patterncycle,'ngait_og1')
    IKdata=importdata('IK_ngait_og1_withfluoro.mot');
elseif strcmp(patterncycle,'ngait_og2')
    IKdata=importdata('IK_ngait_og2_with_fluoro.mot');
elseif strcmp(patterncycle,'ngait_og3')
    IKdata=importdata('IK_ngait_og3_with_fluoro.mot');
elseif strcmp(patterncycle,'2legsquat1')
    IKdata=importdata('IK_2legsquat1_withfluoro.mot');
elseif strcmp(patterncycle,'chairrise1')
    IKdata=importdata('IK_chairrise1_withfluoro.mot');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial1')
    IKdata=importdata('IK_ngait_tm_transition1_trial1_withfluoro.mot');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial2')
    IKdata=importdata('IK_ngait_tm_transition1_trial2_withfluoro.mot');
elseif strcmp(patterncycle,'ngait_tm_transition1_trial3')
    IKdata=importdata('IK_ngait_tm_transition1_trial3_withfluoro.mot');
elseif strcmp(patterncycle,'openfe1')
    IKdata=importdata('IK_openfe1_withfluoro.mot');
elseif strcmp(patterncycle,'bouncy1')
    IKdata=importdata('IK_bouncy1_withfluoro.mot');
elseif strcmp(patterncycle,'medthrust2')
    IKdata=importdata('IK_medthrust2_withfluoro.mot');
elseif strcmp(patterncycle,'moderatecrouch3')
    IKdata=importdata('IK_moderatecrouch3_withfluoro.mot');
else
    keyboard;
end

% Load Moment Arms
if strcmp(patterncycle,'ngait_og4')
    load input_moment_arms_ngait_og4.mat
elseif strcmp(patterncycle,'ngait_og5')
    load input_moment_arms_ngait_og5.mat
elseif strcmp(patterncycle,'ngait_og7')
    load input_moment_arms_ngait_og7.mat
elseif strcmp(patterncycle,'ngait_og1')
    load input_moment_arms_ngait_og1.mat;
elseif strcmp(patterncycle,'ngait_og2')
    load input_moment_arms_ngait_og2.mat;
elseif strcmp(patterncycle,'ngait_og3')
    load input_moment_arms_ngait_og3.mat;
elseif strcmp(patterncycle,'2legsquat1')
    load input_moment_arms_2legsquat1.mat
elseif strcmp(patterncycle,'chairrise1')
    load input_moment_arms_chairrise1.mat;
elseif strcmp(patterncycle,'ngait_tm_transition1_trial1')
    load input_moment_arms_ngait_tm_transition1_trial1.mat
elseif strcmp(patterncycle,'ngait_tm_transition1_trial2')
    load input_moment_arms_ngait_tm_transition1_trial2.mat
elseif strcmp(patterncycle,'ngait_tm_transition1_trial3')
    load input_moment_arms_ngait_tm_transition1_trial3.mat
elseif strcmp(patterncycle,'openfe1')
    load input_moment_arms_openfe1.mat
elseif strcmp(patterncycle,'bouncy1')
    load input_moment_arms_bouncy1.mat
elseif strcmp(patterncycle,'medthrust2')
    load input_moment_arms_medthrust2.mat
elseif strcmp(patterncycle,'moderatecrouch3')
    load input_moment_arms_moderatecrouch3.mat
else
    keyboard;
end
cd(foldermainroot);

[AllID, Allmoment_arms]=Reorder_ID_ma(ID_noCF.ID,moment_arms,Options);

Allmoment_arms_mod=Allmoment_arms+repmat(S.ma_deviation,101,1);

AllID_6toMatch=AllID;
Allmoment_arms_6toMatch=Allmoment_arms_mod;
AllID_6toMatch(:,5:6)=[];
Allmoment_arms_6toMatch(:,(4*44+1):6*44)=[];

%% Load tendon-muscle lengths and velocities
cd(folderExpData);
if strcmp(patterncycle,'ngait_og4')
    load input_lMT_vMT_ngait_og4.mat
elseif strcmp(patterncycle,'ngait_og5')
    load input_lMT_vMT_ngait_og5.mat
elseif strcmp(patterncycle,'ngait_og7')
    load input_lMT_vMT_ngait_og7.mat
elseif strcmp(patterncycle,'ngait_og1')
    load input_lMT_vMT_ngait_og1.mat;
elseif strcmp(patterncycle,'ngait_og2')
    load input_lMT_vMT_ngait_og2.mat
elseif strcmp(patterncycle,'ngait_og3')
    load input_lMT_vMT_ngait_og3.mat
elseif strcmp(patterncycle,'2legsquat1')
    load input_lMT_vMT_2legsquat1.mat
elseif strcmp(patterncycle,'chairrise1')
    load input_lMT_vMT_chairrise1.mat
elseif strcmp(patterncycle,'ngait_tm_transition1_trial1')
    load input_lMT_vMT_ngait_transition1_tm_trial1.mat
elseif strcmp(patterncycle,'ngait_tm_transition1_trial2')
    load input_lMT_vMT_ngait_transition1_tm_trial2.mat
elseif strcmp(patterncycle,'ngait_tm_transition1_trial3')
    load input_lMT_vMT_ngait_transition1_tm_trial3.mat
elseif strcmp(patterncycle,'openfe1')
    load input_lMT_vMT_openfe1.mat
elseif strcmp(patterncycle,'bouncy1')
    load input_lMT_vMT_bouncy1.mat
elseif strcmp(patterncycle,'medthrust2')
    load input_lMT_vMT_medthrust2.mat
elseif strcmp(patterncycle,'moderatecrouch3')
    load input_lMT_vMT_moderatecrouch3.mat
else
    keyboard;
end
cd(foldermainroot);

%% Calculate Rigid Tendon Constants
[c1,c2,lMtilda,FMvtilda] = RigidTendonConstantsVectorized(lMT.data(:,2:(nmuscles+1)),vMT.data(:,2:(nmuscles+1)),params_mod);

%% Calculate Activations
tic
[Activations, IDLoadsMatched, exitflag, Activations_res, Reserve_Actuators]=SolveMuscleActivations(nmuscles,AllID_6toMatch,Allmoment_arms_6toMatch,c1,c2,FM0,Options);
toc

% Solve for Muscle Forces from Activations
[FT,~,~,~,FMpe] = RigidTendonForceVectorized2(Activations,lMT.data(:,2:(nmuscles+1)),vMT.data(:,2:(nmuscles+1)),params_mod);

% Calculate Muscle Forces from c1 and c2 values to verify that they are the same
% MuscleForcesFromCs = Activations.*c1+c2;

for i=1:8
    ID_calculated(:,i)=sum(FT.*(Allmoment_arms(:,((i-1)*44+1):(44*i))),2);
end
% Comparison of medial and lateral contact forces
 [Fymedexp,FyMedFit,Fylatexp,FyLatFit,Txtotexp,Txtotmod]=MedialLateralCalc(moment_arms,Allmoment_arms_mod,FT,ID_noCF,IKdata,patterncycle,Options);