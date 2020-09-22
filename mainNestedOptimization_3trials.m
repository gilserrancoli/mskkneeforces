%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code solves the muscle-force sharing problem using a two-level
% optimization. The outer level optimizes constant parameters. It has the
% option to optimize optimal fiber lengths, tendon slack lengths and muscle
% moment arm offsets. The details of the study can be found in the
% following article:
% Serrancolí, G; Kinney, A; Fregly, B.J. Influence of Musculoskeletal Model 
% Parameter Values on Prediction of Accurate Knee Contact Forces during 
% Walking. Medial Engineering & Physics 2020. In Press. 
% https://doi.org/10.1016/j.medengphy.2020.09.004

% Author: Gil Serrancolí

% The user can change several options at the Options structure. To choose
% optimizing one set of parameters or another, you can change the values of
% Options.opt_lM0, Options.opt_lTs and Options.opt_ma.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mainNestedOptimization_3trials

fprintf('Loading data \n');

%% Specify the name of the walking pattern and the folders
patterncycle1='ngait_og1';
patterncycle2='ngait_og2';
patterncycle3='ngait_og3';
foldermainroot=[pwd '\'];
folderExpData=[foldermainroot, 'ProcessedExperimentalData\'];

%% Specify the Options of the algorithm
Options.match1KneeIDload=0;             %Match 3 Knee ID loads (flexion, adduction, inf-sup) (0), or Match 1 ID load (flexion) (0)
Options.UBactivations=1;                %Whether use muscle activation bounds into Solve for Muscle Activations algorithm (1) or not (0)
Options.trial1=patterncycle1;           %Name of the first trial to analyse
Options.trial2=patterncycle2;           %Name of the second trial to analyse
Options.trial3=patterncycle3;           %Name of the third trial to analyse
Options.StartFromlM0Reas=1;             %Whether start from approx. good initial conditions for lM0 and lTs (1) or literature values (0)
Options.MuscleParam='Arndold';          %Whether use muscle parameters values from Delp 1990 ('Delp') or Arnold 2010 ('Arnold')
Options.expon=5;                        %Exponent for the bound terms at the outer cost function
Options.matchMedialLateral=1;           %Whether to match the Medial Lateral Contact forces (1) or Inferior-Superior force and the Varus-Valgus moment (0).
Options.opt_lM0=1;                      %Whether to optimize lM0 (not equally to lTs) (1) or to set them equal to the literature values)
Options.opt_lTs=1;                      %Whether to optimize lTs (not equally to lM0) (1) or to set them equal to the literature values (0)
Options.opt_ma=1;                       %Whether to optimize moment arm deviations (1) or not (0)
Options.opt_equally_lM0lTs=0;           %Whether to optimize lM0 and lTs with the same scale factor (1) or not (0)
Options.weightingminA_FM0=0;            %Whether use FM0_maxi/FM0 as weighting factors to minimize activations (1) or not (0)
Options.InnerActivationsp05=0;          %Whether minimizing (a_model+0.5)^2 in the Inner Level Cost Function (1) instead of a_model^2 (0)
Options.folderExpData=folderExpData;    %Folder with the experimental data

%% Muscle names
muscles={'addbrev' 'addlong' 'addmagProx' 'addmagMid' 'addmagDist' 'addmagIsch' 'bflh' 'bfsh' 'edl' 'ehl' 'fdl' 'fhl' 'gaslat' 'gasmed' 'gem' 'glmax1' 'glmax2' 'glmax3' 'glmed1' 'glmed2' 'glmed3' 'glmin1' 'glmin2' 'glmin3' 'grac' 'iliacus' 'pect' 'perbrev' 'perlong' 'pertert' 'piri' 'psoas' 'quadfem' 'rf' 'sart' 'semimem' 'semiten' 'soleus' 'tfl' 'tibant' 'tibpost' 'vasint' 'vaslat' 'vasmed'};
nmuscles=length(muscles);

%% Load Inverse Dynamics Loads
[ID1, ID_noCF1, IKdata1, moment_arms1, lMT1, vMT1]=choose_trial(patterncycle1,folderExpData);
[ID2, ID_noCF2, IKdata2, moment_arms2, lMT2, vMT2]=choose_trial(patterncycle2,folderExpData);
[ID3, ID_noCF3, IKdata3, moment_arms3, lMT3, vMT3]=choose_trial(patterncycle3,folderExpData);

if Options.match1KneeIDload
    ID1=ID_noCF1.ID;
    ID2=ID_noCF2.ID;
    ID3=ID_noCF3.ID;
end

[AllID1, Allmoment_arms1]=Reorder_ID_ma(ID1,moment_arms1.moment_arms,Options);
[AllID2, Allmoment_arms2]=Reorder_ID_ma(ID2,moment_arms2.moment_arms,Options);
[AllID3, Allmoment_arms3]=Reorder_ID_ma(ID3,moment_arms3.moment_arms,Options);

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

%% Calculate Normalized Parameter Values for use with the Rigid Tendon Model
params = zeros(34,nmuscles);
[b,c,d,ev,eF,f,g,h] = GetNormalizedParameterValues('exp');
for m = 1:nmuscles
    params(6:40,m) = [b; c; d; ev; eF; f; g; h];
end

params(1,:) = FM0lit;
params(2,:) = lM0lit;
params(3,:) = lTslit;
params(4,:) = PennationAngles*180/pi; % in degrees
params(5,:) = 10*lM0lit; %vMax in m/s (vMax0=10);

%% Calculate the parameters to divide knee sup-inf force into medial and lateral force
[xLat1,xMed1,Fylatexp1,Fymedexp1,Txtot1]=FitEKneeForces(Options.trial1, folderExpData,0);
[xLat2,xMed2,Fylatexp2,Fymedexp2,Txtot2]=FitEKneeForces(Options.trial2, folderExpData,0);
[xLat3,xMed3,Fylatexp3,Fymedexp3,Txtot3]=FitEKneeForces(Options.trial3, folderExpData,0);
Contact_info.xLat1=xLat1;Contact_info.xLat2=xLat2;Contact_info.xLat3=xLat3;
Contact_info.xMed1=xMed1;Contact_info.xMed2=xMed2;Contact_info.xMed3=xMed3;
Contact_info.Fylatexp1=Fylatexp1;Contact_info.Fylatexp2=Fylatexp2;Contact_info.Fylatexp3=Fylatexp3;
Contact_info.Fymedexp1=Fymedexp1;Contact_info.Fymedexp2=Fymedexp2;Contact_info.Fymedexp3=Fymedexp3;
Contact_info.Txtot1=Txtot1;Contact_info.Txtot2=Txtot2;Contact_info.Txtot3=Txtot3;

  
%% Solve for lMT, lTs, ma and FM0, minimizing Passive Muscle Force and Muscle Activations and taking into account Activation Boundaries
params_mod=params;
params_mod(1,:)=2*params(1,:);

aux=load('input_lM0lTs_Optimized_to_maxlMtilda1_3gaitTrials.mat');
params_mod(2,:)=lM0lit.*aux.slM0;
params_mod(3,:)=lTslit.*aux.slTs;

pos_nonzero_ma=find(sum(Allmoment_arms1)~=0);
if Options.match1KneeIDload
    ma_deviation_x=zeros(1,size(AllID1,2)*44-44*2);
    pos_nonzero_ma((pos_nonzero_ma>44*4)&(pos_nonzero_ma<=44*6))=[];
    pos_nonzero_ma(pos_nonzero_ma>44*6)=pos_nonzero_ma(pos_nonzero_ma>44*6)-44*2;
else
    ma_deviation_x=zeros(1,size(AllID1,2)*44);
end

x0=[];lb=[];ub=[];
if Options.opt_lM0
    if Options.StartFromlM0Reas
        x0(1:44)=aux.slM0;
    else
        x0(1:44)=1;
    end
     lb=zeros(1,nmuscles);ub=3*ones(1,nmuscles);
end
if Options.opt_lTs
    if Options.StartFromlM0Reas
        x0=[x0 aux.slM0];
    else
        x0=[x0 ones(1,44)];
    end
     lb=[lb zeros(1,44)];ub=[ub 3*ones(1,44)];
end
if (isempty(x0)&&Options.opt_equally_lM0lTs)
    if Options.StartFromlM0Reas
        x0=mean([aux.slM0; aux.slTs]);
    else
        x0=ones(1,44);
    end
    lb=zeros(1,nmuscles);ub=3*ones(1,nmuscles);
end    
if Options.opt_ma
      x0=[x0, ma_deviation_x(pos_nonzero_ma)]; %scale factor ma
      lb=[lb -20*ones(1,147)];ub=[ub 20*ones(1,147)];
end
    
p_before_ty=find(pos_nonzero_ma<=5*44);
p_ty=find((pos_nonzero_ma>5*44)&(pos_nonzero_ma<=6*44));
p_after_ty=find(pos_nonzero_ma>6*44);
if Options.opt_lTs && ~Options.opt_lM0 && ~Options.opt_ma
    lb(1:44)=1e-1;
elseif Options.opt_lTs && Options.opt_lM0 && Options.opt_ma
    lb(1:88)=1e-1;
    lb((88+1):(88+length(pos_nonzero_ma)))=-1.2;
    ub((88+1):(88+length(pos_nonzero_ma)))=1.2;

end

options_lsqnonlin=optimset('Display','iter','UseParallel','always','MaxFunEvals',100000,'OutputFcn',@SaveIter);tic
fprintf('Running optimization \n');
t0=tic;
[x,resnorm,residual,exitflag,output]=lsqnonlin(@costfunOuterLevel_3gaittrials,x0,lb,ub,options_lsqnonlin,params_mod,params,lMT1,lMT2,lMT3,vMT1,vMT2,vMT3,nmuscles,Allmoment_arms1,Allmoment_arms2,Allmoment_arms3,AllID1,AllID2,AllID3,pos_nonzero_ma,muscles,Contact_info,ID_noCF1,ID_noCF2,ID_noCF3,Options);
toc(t0)

    
%% Recalculate lM0 and lTs
if Options.match1KneeIDload
    ma_deviation=zeros(1,size(AllID1,2)*44-44*2);
    ma_deviation_x=zeros(1,size(AllID1,2)*44-44*2);    
else
    ma_deviation=zeros(1,size(AllID1,2)*44);
    ma_deviation_x=zeros(1,size(AllID1,2)*44);
end
            
if Options.opt_lM0
    lM0=lM0lit.*x(1:44);
    params_mod(2,:)=lM0;
    if Options.opt_lTs
        lTs=lTslit.*x(45:88);
        params_mod(3,:)=lTs;
        if Options.opt_ma    
            ma_deviation_x(pos_nonzero_ma)=x((88+1):(88+length(pos_nonzero_ma)));
        else
        end
    else
        lTs=lTslit;
        if Options.opt_ma
            %Calculate current moment arms
            ma_deviation_x(pos_nonzero_ma)=x((44+1):(44+length(pos_nonzero_ma)));
        else
        end
    end
else
    lM0=lM0lit;
    if Options.opt_lTs
        lTs=lTslit.*x(1:44);
        params_mod(3,:)=lTs;
        if Options.opt_ma    
            ma_deviation_x(pos_nonzero_ma)=x((44+1):(44+length(pos_nonzero_ma)));
        else
        end
    else
        lTs=lTslit;
        if Options.opt_equally_lM0lTs
            lM0=lM0lit.*x(1:44);
            lTs=lTslit.*x(1:44);
            if Options.opt_ma
                %Calculate current moment arms
                ma_deviation_x(pos_nonzero_ma)=x(45:(44+length(pos_nonzero_ma)));
            else
            end
        else
            if Options.opt_ma
                %Calculate current moment arms
                ma_deviation_x(pos_nonzero_ma)=x(1:length(pos_nonzero_ma));
            else
            end
        end
    end
end

if Options.match1KneeIDload
    ma_deviation=0.005*ma_deviation_x;
else
    p_before_ty=find(pos_nonzero_ma<=5*44);
    p_ty=find((pos_nonzero_ma>5*44)&(pos_nonzero_ma<=6*44));
    p_after_ty=find(pos_nonzero_ma>6*44);
    ma_deviation(:,pos_nonzero_ma([p_before_ty p_after_ty]))=0.005*ma_deviation_x(:,pos_nonzero_ma([p_before_ty p_after_ty]));
    ma_deviation(:,pos_nonzero_ma(p_ty))=0.015*ma_deviation_x(:,pos_nonzero_ma(p_ty));
end
if Options.match1KneeIDload
    Allmoment_arms_mod1(:,1:44*4)=Allmoment_arms1(:,1:44*4)+repmat(ma_deviation(1:44*4),101,1);
    Allmoment_arms_mod1(:,(44*4+1):(44*6))=Allmoment_arms1(:,(44*4+1):(44*6));
    Allmoment_arms_mod1(:,(44*6+1):(44*8))=Allmoment_arms1(:,(44*6+1):(44*8))+repmat(ma_deviation((44*4+1):44*6),101,1);
    Allmoment_arms_mod2(:,1:44*4)=Allmoment_arms2(:,1:44*4)+repmat(ma_deviation(1:44*4),101,1);
    Allmoment_arms_mod2(:,(44*4+1):(44*6))=Allmoment_arms2(:,(44*4+1):(44*6));
    Allmoment_arms_mod2(:,(44*6+1):(44*8))=Allmoment_arms2(:,(44*6+1):(44*8))+repmat(ma_deviation((44*4+1):44*6),101,1);
    Allmoment_arms_mod3(:,1:44*4)=Allmoment_arms3(:,1:44*4)+repmat(ma_deviation(1:44*4),101,1);
    Allmoment_arms_mod3(:,(44*4+1):(44*6))=Allmoment_arms3(:,(44*4+1):(44*6));
    Allmoment_arms_mod3(:,(44*6+1):(44*8))=Allmoment_arms3(:,(44*6+1):(44*8))+repmat(ma_deviation((44*4+1):44*6),101,1);
else
    Allmoment_arms_mod1=Allmoment_arms1+repmat(ma_deviation,101,1);
    Allmoment_arms_mod2=Allmoment_arms2+repmat(ma_deviation,101,1);
    Allmoment_arms_mod3=Allmoment_arms3+repmat(ma_deviation,101,1);
end

%% Recalculate FM0
FM0=params_mod(1,:);

%% Calculate c1 and c2, and modified moment arms
[c1_og1,c2_og1,lMtilda1,FMvtilda1] = RigidTendonConstantsVectorized(lMT1.data(:,2:(nmuscles+1)),vMT1.data(:,2:(nmuscles+1)),params_mod);
[c1_og2,c2_og2,lMtilda2,FMvtilda2] = RigidTendonConstantsVectorized(lMT2.data(:,2:(nmuscles+1)),vMT2.data(:,2:(nmuscles+1)),params_mod);
[c1_og3,c2_og3,lMtilda3,FMvtilda3] = RigidTendonConstantsVectorized(lMT3.data(:,2:(nmuscles+1)),vMT3.data(:,2:(nmuscles+1)),params_mod);

AllID_6toMatch1=AllID1;
AllID_6toMatch2=AllID2;
AllID_6toMatch3=AllID3;
Allmoment_arms_mod_6toMatch1=Allmoment_arms_mod1;
Allmoment_arms_mod_6toMatch2=Allmoment_arms_mod2;    
Allmoment_arms_mod_6toMatch3=Allmoment_arms_mod3;    
AllID_6toMatch1(:,5:6)=[];
AllID_6toMatch2(:,5:6)=[];
AllID_6toMatch3(:,5:6)=[];
Allmoment_arms_mod_6toMatch1(:,[(4*44+1):6*44])=[];
Allmoment_arms_mod_6toMatch2(:,[(4*44+1):6*44])=[];
Allmoment_arms_mod_6toMatch3(:,[(4*44+1):6*44])=[];

%% Solve for Muscle Activations 
[Activations1, IDLoadsMatched1, exitflag1, Activations_res1, Reserve_Actuators1]=SolveMuscleActivations(nmuscles,AllID_6toMatch1,Allmoment_arms_mod_6toMatch1,c1_og1,c2_og1,FM0,Options);
[Activations2, IDLoadsMatched2, exitflag2, Activations_res2, Reserve_Actuators2]=SolveMuscleActivations(nmuscles,AllID_6toMatch2,Allmoment_arms_mod_6toMatch2,c1_og2,c2_og2,FM0,Options);
[Activations3, IDLoadsMatched3, exitflag3, Activations_res3, Reserve_Actuators3]=SolveMuscleActivations(nmuscles,AllID_6toMatch3,Allmoment_arms_mod_6toMatch3,c1_og3,c2_og3,FM0,Options);

% Solve for Muscle Forces from Activations
[FT1,~,~,~,FMpe1] = RigidTendonForceVectorized2(Activations1,lMT1.data(:,2:(nmuscles+1)),vMT1.data(:,2:(nmuscles+1)),params_mod);
[FT2,~,~,~,FMpe2] = RigidTendonForceVectorized2(Activations2,lMT2.data(:,2:(nmuscles+1)),vMT2.data(:,2:(nmuscles+1)),params_mod);
[FT3,~,~,~,FMpe3] = RigidTendonForceVectorized2(Activations3,lMT3.data(:,2:(nmuscles+1)),vMT3.data(:,2:(nmuscles+1)),params_mod);

% Calculate Muscle Forces from c1 and c2 values to verify that they are the same
MuscleForcesFromCs1 = Activations1.*c1_og1+c2_og1;
MuscleForcesFromCs2 = Activations2.*c1_og2+c2_og2;
MuscleForcesFromCs3 = Activations3.*c1_og3+c2_og3;

% Comparison of ID load calculated and experimental ID loads
if Options.match1KneeIDload
    nloads_toMatch=6;
else
    nloads_toMatch=8;
end
for i=1:nloads_toMatch
    ID_calculated1(:,i)=sum(FT1.*(Allmoment_arms1(:,((i-1)*44+1):(44*i))),2);
    ID_calculated2(:,i)=sum(FT2.*(Allmoment_arms2(:,((i-1)*44+1):(44*i))),2);
    ID_calculated3(:,i)=sum(FT3.*(Allmoment_arms3(:,((i-1)*44+1):(44*i))),2);
end

% Comparison of medial and lateral contact forces
[Fymedexp1,FyMedFit1,Fylatexp1,FyLatFit1,Txtotexp1,Txtotmod1]=MedialLateralCalc(moment_arms1.moment_arms,Allmoment_arms_mod1,FT1,ID_noCF1,IKdata1,patterncycle1,Options);
[Fymedexp2,FyMedFit2,Fylatexp2,FyLatFit2,Txtotexp2,Txtotmod2]=MedialLateralCalc(moment_arms2.moment_arms,Allmoment_arms_mod2,FT2,ID_noCF2,IKdata2,patterncycle2,Options);
[Fymedexp3,FyMedFit3,Fylatexp3,FyLatFit3,Txtotexp3,Txtotmod3]=MedialLateralCalc(moment_arms3.moment_arms,Allmoment_arms_mod3,FT3,ID_noCF3,IKdata3,patterncycle3,Options);

keyboard;
end

function [ID, ID_noCF, IKdata, moment_arms, lMT, vMT]=choose_trial(pattern,expdata_folder)
current_folder=pwd;
cd(expdata_folder)
switch    pattern
    case 'ngait_og1'
        ID_wCF=load('input_ID_withContactForces_ngait_og1.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og1.mat');  
        IKdata=importdata('IK_ngait_og1_withfluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og1.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og1.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_og2'
        ID_wCF=load('input_ID_withContactForces_ngait_og2.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og2.mat');
        IKdata=importdata('IK_ngait_og2_with_fluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og2.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og2.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_og3'  
        ID_wCF=load('input_ID_withContactForces_ngait_og3.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og3.mat'); 
        IKdata=importdata('IK_ngait_og3_with_fluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og3.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og3.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_og4'
        ID_wCF=load('input_ID_withContactForces_ngait_og4.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og4.mat');
        IKdata=importdata('IK_ngait_og4_with_fluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og4.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og4.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_og5'
        ID_wCF=load('input_ID_withContactForces_ngait_og5.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og5.mat');
        IKdata=importdata('IK_ngait_og5_Extmuscle_withfluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og5.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og5.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_og7'
        ID_wCF=load('input_ID_withContactForces_ngait_og7.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_og7.mat');
        IKdata=importdata('IK_ngait_og7_with_fluoro.mot');
        moment_arms=load('input_moment_arms_ngait_og7.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_og7.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_tm_transition1_trial1'
        ID_wCF=load('input_ID_withContactForces_ngait_tm_transition1_trial1.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial1.mat');
        IKdata=importdata('IK_ngait_tm_transition1_trial1_withfluoro.mot');
        moment_arms=load('input_moment_arms_ngait_tm_transition1_trial1.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_transition1_tm_trial1.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_tm_transition1_trial2'
        ID_wCF=load('input_ID_withContactForces_ngait_tm_transition1_trial2.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial2.mat');
        IKdata=importdata('IK_ngait_tm_transition1_trial2_withfluoro.mot');
        moment_arms=load('input_moment_arms_ngait_tm_transition1_trial2.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_transition1_tm_trial2.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'ngait_tm_transition1_trial3'
        ID_wCF=load('input_ID_withContactForces_ngait_tm_transition1_trial3.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_ngait_tm_transition1_trial3.mat');
        IKdata=importdata('IK_ngait_tm_transition1_trial3_withfluoro.mot');
        moment_arms=load('input_moment_arms_ngait_tm_transition1_trial3.mat');
        lMTvMT1=load('input_lMT_vMT_ngait_transition1_tm_trial3.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'chairrise1'
        ID_wCF=load('input_ID_withContactForces_chairrise1.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_chairrise1.mat');
        IKdata=importdata('IK_chairrise1_withfluoro.mot');
        moment_arms=load('input_moment_arms_chairrise1.mat');
        lMTvMT1=load('input_lMT_vMT_chairrise1.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case '2legsquat1'
        ID_wCF=load('input_ID_withContactForces_2legsquat1.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_2legsquat1.mat');  
        IKdata=importdata('IK_2legsquat1_withfluoro.mot');  
        moment_arms=load('input_moment_arms_2legsquat1.mat');
        lMTvMT1=load('input_lMT_vMT_2legsquat1.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
    case 'bouncy1'
        ID_wCF=load('input_ID_withContactForces_bouncy1.mat');ID=ID_wCF.ID;
        ID_noCF=load('input_ID_withoutContactForces_bouncy1.mat');  
        IKdata=importdata('IK_bouncy1_withfluoro.mot');  
        moment_arms=load('input_moment_arms_bouncy1.mat');
        lMTvMT1=load('input_lMT_vMT_bouncy1.mat');lMT=lMTvMT1.lMT;vMT=lMTvMT1.vMT;
end
cd(current_folder);
end