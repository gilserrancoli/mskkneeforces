function [FyLatFit,FyMedFit]=Calculate_MedialLateralFit(moment_arms_mod,FM,ID_noCF,xLat,xMed)

% Author: Gil Serrancolí

%% Calculate the medial and lateral contact forces from the inf-sup contact force and the varus-valgus contact moment

% inf. sup force
ma=moment_arms_mod(:,(5*44+1):6*44);
MuscleIDloads(:,1)=sum(FM.*ma,2); 
%knee adduction
ma=moment_arms_mod(:,(4*44+1):5*44);
MuscleIDloads(:,2)=sum(FM.*ma,2);   

AllIDloads=ID_noCF.ID.data; 
TotalIDloads(:,1)=AllIDloads(:,15);%F knee_ty 
TotalIDloads(:,2)=AllIDloads(:,11);%M knee_adduction
CFIDloads_knee_onfemur=TotalIDloads-MuscleIDloads;
CFIDloads_knee_ontibia=-CFIDloads_knee_onfemur;
CFIDloads_tibialtray=CFIDloads_knee_ontibia;
A=[CFIDloads_tibialtray(:,1) CFIDloads_tibialtray(:,2)]; %knee_ty and knee_adduction
A(:,2)=A(:,2)*1000;
FyLatFit=A*xLat;
FyMedFit=A*xMed;


end