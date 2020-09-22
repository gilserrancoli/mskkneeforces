function [xLat,xMed,FyLat,FyMed,TxTot]=FitEKneeForces(gait_trial, folder_ExpData, bplot)

close all

% Read in contact force data from pose estimation optimization
% Data are stored in the following columns:
% Column 1: Frame
% Column 2: FxMed
% Column 3: FyMed
% Column 4: FzMed
% Column 5: TxMed
% Column 6: TyMed
% Column 7: TzMed
% Column 8: FxLat
% Column 9: FyLat
% Column 10: FzLat
% Column 11: TxLat
% Column 12: TyLat
% Column 13: TzLat
% Column 14: FxTot
% Column 15: FyTot
% Column 16: FzTot
% Column 17: TxTot
% Column 18: TyTot
% Column 19: TzTot

current_folder=pwd;
cd(folder_ExpData);

if strcmp(gait_trial,'ngait_og1')
    Data=xlsread('Contact Loads for JW ngait og1.xlsx','A3:S103');
elseif strcmp(gait_trial,'ngait_og2')
    Data=xlsread('Contact Loads for JW ngait og2.xlsx','A3:S103');
elseif strcmp(gait_trial,'ngait_og3')
    Data=xlsread('Contact Loads for JW ngait og3.xlsx','A3:S103');
elseif strcmp(gait_trial,'ngait_og4')
    Data=xlsread('Contact Loads for JW ngait og4.xlsx','A3:S103');
elseif strcmp(gait_trial,'ngait_og5')
    Data=xlsread('Contact Loads for JW ngait og5.xlsx','A3:S103');
elseif strcmp(gait_trial,'ngait_og7')
    Data=xlsread('Contact Loads for JW ngait og7.xlsx');
elseif strcmp(gait_trial,'ngait_tm_transition1_trial1')
    impdata=importdata('Contact Loads for JW ngait tm transition1_trial1_smooth.xlsx');
    Data=impdata.data.Data;
elseif strcmp(gait_trial,'ngait_tm_transition1_trial2')
    impdata=importdata('Contact Loads for JW ngait tm transition1_trial2_smooth.xlsx');
    Data=impdata.data.Data;
elseif strcmp(gait_trial,'ngait_tm_transition1_trial3')
    impdata=importdata('Contact Loads for JW ngait tm transition1_trial3_smooth.xlsx');
    Data=impdata.data.Data;
elseif strcmp(gait_trial,'chairrise1')
    Data=xlsread('Contact Loads for JW chairrise1_smooth.xlsx');
elseif strcmp(gait_trial,'2legsquat1')
    Data=xlsread('Contact Loads for JW 2legsquat1_smooth.xlsx');
elseif strcmp(gait_trial,'openfe1')
    Data=xlsread('Contact Loads for JW openfe1.xlsx');
elseif strcmp(gait_trial,'bouncy1')
    impdata=importdata('Contact Loads for JW bouncy1_smooth.xlsx');
    Data=impdata.data.Data;
elseif strcmp(gait_trial,'medthrust2')
    impdata=importdata('Contact Loads for JW medthrust2_smooth.xlsx');
    Data=impdata.data.Data;
elseif strcmp(gait_trial,'moderatecrouch3')
    impdata=importdata('Contact Loads for JW moderatecrouch3_smooth.xlsx');
    Data=impdata.data.Data;
else
    keyboard;
end

cd(current_folder);

npts = size(Data,1);
Percent = linspace(0,100,npts)';

% Extract contact force data for regression fitting
FyMed = Data(:,3);
FyLat = Data(:,9);
FyTot = Data(:,15);
TxTot = Data(:,17);
TzTot = Data(:,19);

% Fit medial contact force as a function of FyTot and TxTot
% FyMed = x(1,1)*FyTot+x(2,1)*TxTot
A = [FyTot TxTot];
b = FyMed;
xMed = A\b;

FyMedFit = A*xMed;

if bplot
    subplot(1,2,1), plot(Percent,FyMed,'ko','LineWidth',2)
    hold on
    subplot(1,2,1), plot(Percent,FyMedFit,'r-','LineWidth',2)
    legend('Experiment','Regression')
    xlabel('Gait Cycle (%)')
    ylabel('Medial Contact Force (N)')
end

% Fit lateral contact force as a function of FyTot and TxTot
% FyLat = x(1,1)*FyTot+x(2,1)*TxTot
A = [FyTot TxTot];
b = FyLat;
xLat = A\b;

FyLatFit = A*xLat;

if bplot
    subplot(1,2,2), plot(Percent,FyLat,'ko','LineWidth',2)
    hold on
    subplot(1,2,2), plot(Percent,FyLatFit,'r-','LineWidth',2)
    legend('Experiment','Regression')
    xlabel('Gait Cycle (%)')
    ylabel('Lateral Contact Force (N)')
end
