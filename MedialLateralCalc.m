function [Fymedexp,FyMedFit,Fylatexp,FyLatFit,Txtot_exp,Txtot_mod]=MedialLateralCalc(moment_arms,moment_arms_mod,FM,ID_noCF,IKdata,gait_trial,Options)

% Author: Gil Serrancolí

%% Calculate medial and lateral contact forces from muscle moments and inverse dynamics loads

for i=1:size(moment_arms,2);
        loadID=moment_arms(i).coordinate;
        if strcmp(moment_arms(i).coordinate,loadID)
            ma=moment_arms(i).data(:,2:end);
            if strcmp(loadID,'knee_tx')
                MuscleIDloads(:,1)=sum(FM.*ma,2);
            elseif strcmp(loadID,'knee_ty')
                if Options.match1KneeIDload
                else
                    ma=moment_arms_mod(:,(5*44+1):6*44);
                end
                    MuscleIDloads(:,2)=sum(FM.*ma,2);
            elseif strcmp(loadID,'knee_tz')
                MuscleIDloads(:,3)=sum(FM.*ma,2);
            elseif strcmp(loadID,'knee_adduction')
                if Options.match1KneeIDload
                else
                    ma=moment_arms_mod(:,(4*44+1):5*44);
                end
                MuscleIDloads(:,4)=sum(FM.*ma,2);    
            elseif strcmp(loadID,'knee_rotation')
                MuscleIDloads(:,5)=sum(FM.*ma,2);   
            elseif strcmp(loadID,'knee_flexion')
                ma=moment_arms_mod(:,(3*44+1):4*44);
                MuscleIDloads(:,6)=sum(FM.*ma,2);   
            else
            end
        end
end

AllIDloads=ID_noCF.ID.data;
TotalIDloads(:,1)=AllIDloads(:,14);%F knee_tx
TotalIDloads(:,2)=AllIDloads(:,15);%F knee_ty
TotalIDloads(:,3)=AllIDloads(:,16);%F knee_tz
TotalIDloads(:,4)=AllIDloads(:,11);%M knee_adduction
TotalIDloads(:,5)=AllIDloads(:,12);%M knee_rootation
TotalIDloads(:,6)=AllIDloads(:,13);%M knee_flexion

CFIDloads_knee_onfemur=TotalIDloads-MuscleIDloads;
CFIDloads_knee_ontibia=-CFIDloads_knee_onfemur;

% translations_raw(:,1)=IKdata.data(:,14); %knee_tx
% translations_raw(:,2)=IKdata.data(:,15); %knee_ty
% translations_raw(:,3)=IKdata.data(:,16); %knee_tz

CFIDloads_tibialtray(:,1:3)=CFIDloads_knee_ontibia(:,1:3);
CFIDloads_tibialtray(:,4:6)=CFIDloads_knee_ontibia(:,4:6);%+cross(translations,CFIDloads_knee_ontibia(:,1:3));

figure(1);
[xLat,xMed,Fylatexp,Fymedexp,Txtot]=FitEKneeForces(gait_trial,Options.folderExpData,1);

A=[CFIDloads_tibialtray(:,2) CFIDloads_tibialtray(:,4)]; %knee_ty and knee_adduction
A(:,2)=A(:,2)*1000;
FyLatFit=A*xLat;
FyMedFit=A*xMed;

figure(2);
npts=length(FyLatFit);
Percent = linspace(0,100,npts)';
subplot(1,2,1), plot(Percent,-FyMedFit,'k','LineWidth',2);hold all
plot(Percent,-Fymedexp,'r','LineWidth',2);
legend('Model','Experimental')
xlabel('Gait Cycle (%)')
R2=1-(sum((Fymedexp-FyMedFit).^2)./(sum((Fymedexp-mean(Fymedexp)).^2)));
RMSE=sqrt(sum((Fymedexp-FyMedFit).^2)/length(Fymedexp));
title(['Medial Contact Force (N)',10,'R^2=',num2str(R2),10,'RMSE=',num2str(RMSE)])

subplot(1,2,2), plot(Percent,-FyLatFit,'k','LineWidth',2);hold all
plot(Percent,-Fylatexp,'r','LineWidth',2);
legend('Model','Experimental')
xlabel('Gait Cycle (%)');
R2=1-(sum((Fylatexp-FyLatFit).^2)./(sum((Fylatexp-mean(Fylatexp)).^2)));
RMSE=sqrt(sum((Fylatexp-FyLatFit).^2)/length(Fylatexp));
title(['Lateral Contact Force (N)',10,'R^2=',num2str(R2),10,'RMSE=',num2str(RMSE)]);


figure(3)
plot(Percent,CFIDloads_tibialtray(:,4),'k','LineWidth',2);hold all
plot(Percent,Txtot/1000,'r','LineWidth',2);
legend('Model','Experimental');
xlabel('Gait Cycle (%)');
R2=1-(sum((Txtot/1000-CFIDloads_tibialtray(:,4)).^2)./(sum((Txtot/1000-mean(Txtot/1000)).^2)));  
RMSE=sqrt(sum((Txtot/1000-CFIDloads_tibialtray(:,4)).^2)/length(Txtot));
title(['Knee Adduction Contact Moment (Nm)',10,'R^2=',num2str(R2),10,'RMSE=',num2str(RMSE)]);
Txtot_exp=Txtot;
Txtot_mod=CFIDloads_tibialtray(:,4);

figure(4);
plot(Percent,-FyMedFit-FyLatFit,'k','LineWidth',2);hold all;
plot(Percent,-Fymedexp-Fylatexp,'r','LineWidth',2);
legend('Model','Experimental');
xlabel('Gait Cycle (%)');
R2=1-(sum((Fylatexp+Fymedexp-FyLatFit-FyMedFit).^2)./(sum((Fylatexp+Fymedexp-mean(Fylatexp+Fymedexp)).^2)));
RMSE=sqrt(sum((Fylatexp+Fymedexp-FyLatFit-FyMedFit).^2)/length(Fylatexp+Fymedexp));
title(['Inf.-Sup. Contact force (N)',10,'R^2=',num2str(R2),10,'RMSE=',num2str(RMSE)]);
end