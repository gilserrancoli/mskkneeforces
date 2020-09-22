function f=costfunOuterLevel_3gaittrials(x,params_mod,paramslit,lMT1,lMT2,lMT3,vMT1,vMT2,vMT3,nmuscles,Allmoment_arms1,Allmoment_arms2,Allmoment_arms3,AllID1,AllID2,AllID3,pos_nonzero_ma,muscles,Contact_info,ID_noCF1,ID_noCF2,ID_noCF3,Options)
% Author: Gil Serrancolí

%% cost function of the inner-level optimization

lM0lit=paramslit(2,:);
lTslit=paramslit(3,:);
expon=Options.expon;

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
        if Options.opt_ma
            %Calculate current moment arms
            ma_deviation_x(pos_nonzero_ma)=x(1:length(pos_nonzero_ma));
        else
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
    
%% Calculate current FM0
FM0=params_mod(1,:);

[c1_og1,c2_og1,lMtilda1,FMvtilda1] = RigidTendonConstantsVectorized(lMT1.data(:,2:(nmuscles+1)),vMT1.data(:,2:(nmuscles+1)),params_mod);
[c1_og2,c2_og2,lMtilda2,FMvtilda2] = RigidTendonConstantsVectorized(lMT2.data(:,2:(nmuscles+1)),vMT2.data(:,2:(nmuscles+1)),params_mod);
[c1_og3,c2_og3,lMtilda3,FMvtilda3] = RigidTendonConstantsVectorized(lMT3.data(:,2:(nmuscles+1)),vMT3.data(:,2:(nmuscles+1)),params_mod);


%% Solve for Muscle Activations
AllID_6toMatch1=AllID1;
AllID_6toMatch2=AllID2;
AllID_6toMatch3=AllID3;
Allmoment_arms_mod_6toMatch1=Allmoment_arms_mod1;
Allmoment_arms_mod_6toMatch2=Allmoment_arms_mod2;
Allmoment_arms_mod_6toMatch3=Allmoment_arms_mod3;
AllID_6toMatch1(:,5:6)=[];
AllID_6toMatch2(:,5:6)=[];
AllID_6toMatch3(:,5:6)=[];
Allmoment_arms_mod_6toMatch1(:,(4*44+1):6*44)=[];
Allmoment_arms_mod_6toMatch2(:,(4*44+1):6*44)=[];
Allmoment_arms_mod_6toMatch3(:,(4*44+1):6*44)=[];


try

% Activations as design variables
[Activations1, IDLoadsMatched1, exitflag1, Activations_res1, Reserve_Actuators1]=SolveMuscleActivations(nmuscles,AllID_6toMatch1,Allmoment_arms_mod_6toMatch1,c1_og1,c2_og1,FM0,Options);
[Activations2, IDLoadsMatched2, exitflag2, Activations_res2, Reserve_Actuators2]=SolveMuscleActivations(nmuscles,AllID_6toMatch2,Allmoment_arms_mod_6toMatch2,c1_og2,c2_og2,FM0,Options);
[Activations3, IDLoadsMatched3, exitflag3, Activations_res3, Reserve_Actuators3]=SolveMuscleActivations(nmuscles,AllID_6toMatch3,Allmoment_arms_mod_6toMatch3,c1_og3,c2_og3,FM0,Options);


catch
    keyboard;
end

[FT1,~,~,~,~] = RigidTendonForceVectorized2(Activations1,lMT1.data(:,2:(nmuscles+1)),vMT1.data(:,2:(nmuscles+1)),params_mod);
[FT2,~,~,~,~] = RigidTendonForceVectorized2(Activations2,lMT2.data(:,2:(nmuscles+1)),vMT2.data(:,2:(nmuscles+1)),params_mod);
[FT3,~,~,~,~] = RigidTendonForceVectorized2(Activations3,lMT3.data(:,2:(nmuscles+1)),vMT3.data(:,2:(nmuscles+1)),params_mod);

f=[];

%% Minimization of reserve activations

f=[f; reshape((Activations_res1.^2),101*6,1)];
f=[f; reshape((Activations_res2.^2),101*6,1)];
f=[f; reshape((Activations_res3.^2),101*6,1)];


%% Scale factors deviation

if Options.opt_lM0||Options.opt_lTs
    sclM0=lM0./lM0lit;
    sclTs=lTs./lTslit;
    f=[f; (((sclM0-sclTs)./(0.2*mean([sclM0;sclTs])))).^expon'];
else
end
if Options.opt_ma
     f=[f; (ma_deviation(:,pos_nonzero_ma(p_ty))'/0.015).^expon];
     f=[f; (ma_deviation(:,pos_nonzero_ma([p_before_ty p_after_ty]))'/0.005).^expon];
end

%% lMtilda close to 1
k=500;
smooth_max1=log(sum(exp(k*lMtilda1)))/k;
aux_smooth1=1/2+(1/2)*tanh(100*(smooth_max1-1));
smooth_max2=log(sum(exp(k*lMtilda2)))/k;
aux_smooth2=1/2+(1/2)*tanh(100*(smooth_max2-1));
smooth_max3=log(sum(exp(k*lMtilda3)))/k;
aux_smooth3=1/2+(1/2)*tanh(100*(smooth_max3-1));
f=[f; aux_smooth1'.*reshape(((smooth_max1-1)/0.2).^expon,44,1)];
f=[f; aux_smooth2'.*reshape(((smooth_max2-1)/0.2).^expon,44,1)];
f=[f; aux_smooth3'.*reshape(((smooth_max3-1)/0.2).^expon,44,1)];


%% Tracking errors 

if Options.match1KneeIDload
    
else
   % Tracking medial and lateral loads
   if Options.matchMedialLateral
        [FyLatFit1,FyMedFit1]=Calculate_MedialLateralFit(Allmoment_arms_mod1,FT1,ID_noCF1,Contact_info.xLat1,Contact_info.xMed1);
        [FyLatFit2,FyMedFit2]=Calculate_MedialLateralFit(Allmoment_arms_mod2,FT2,ID_noCF2,Contact_info.xLat2,Contact_info.xMed2);
        [FyLatFit3,FyMedFit3]=Calculate_MedialLateralFit(Allmoment_arms_mod3,FT3,ID_noCF3,Contact_info.xLat3,Contact_info.xMed3);
        f=[f; ((Contact_info.Fylatexp1-FyLatFit1)/20)];
        f=[f; ((Contact_info.Fylatexp2-FyLatFit2)/20)];
        f=[f; ((Contact_info.Fylatexp3-FyLatFit3)/20)];
        
        f=[f; ((Contact_info.Fymedexp1-FyMedFit1)/20)];
        f=[f; ((Contact_info.Fymedexp2-FyMedFit2)/20)];
        f=[f; ((Contact_info.Fymedexp3-FyMedFit3)/20)];
   else
       % Tracking knee inf-sup. forces and knee adduction moment
        f=[f; (AllID1(:,5)-sum(FT1.*Allmoment_arms_mod1(:,(4*44+1):5*44),2))/0.5];
        f=[f; (AllID2(:,5)-sum(FT2.*Allmoment_arms_mod2(:,(4*44+1):5*44),2))/0.5];
        f=[f; (AllID3(:,5)-sum(FT3.*Allmoment_arms_mod3(:,(4*44+1):5*44),2))/0.5];

        f=[f; (AllID1(:,6)-sum(FT1.*Allmoment_arms_mod1(:,(5*44+1):6*44),2))/20];
        f=[f; (AllID2(:,6)-sum(FT2.*Allmoment_arms_mod2(:,(5*44+1):6*44),2))/20];
        f=[f; (AllID3(:,6)-sum(FT3.*Allmoment_arms_mod3(:,(5*44+1):6*44),2))/20];
   end
end

end


