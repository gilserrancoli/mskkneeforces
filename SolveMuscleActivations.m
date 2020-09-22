function [Activations, IDLoadsMatched, exitflag, Activations_res, Reserve_Actuators]=SolveMuscleActivations(nmuscles,AllID,Allmoment_arms,c1,c2,FM0,Options)
% Author: Gil Serrancolí

%% Solve for muscle activations with quadprog

Activations_ext=[];
exitflag=[];
IDLoadsMatched=[];
H = 2*(1)*eye(nmuscles+6,nmuscles+6);
if Options.weightingminA_FM0
   H(1:nmuscles,1:nmuscles)=diag(FM0./sum(FM0)); 
end

f = zeros(nmuscles+6,1);
if Options.InnerActivationsp05
    f(1:nmuscles)=1;
end
lb = 1e-5*ones(nmuscles,1);
lb=[lb;  -ones(6,1)*inf];
if Options.UBactivations
    ub=ones(nmuscles,1);
    ub=[ub; ones(6,1)*inf];
else
    ub = [];
end
T=0.5;
opts = optimset('Display','off');
warning('off');

parfor j=1:11
Activations_ext_aux=[];
IDLoadsMatched_aux=[];
exitflag_aux=[];
for k=1:5
    t=2*((j-1)*5+k)-1;
    if t>101
    else
        % Form Aeq matrix of c1*moment arms (8xnmuscles)
        Aeq=zeros(size(AllID,2),44+6);
        sum_c2=zeros(size(AllID,2),44);
        for m = 1:nmuscles
            Aeq(:,m) = (c1(t,m)*Allmoment_arms(t,m:nmuscles:end))';
            sum_c2(:,m) = (c2(t,m)*Allmoment_arms(t,m:nmuscles:end))';
        end

        Aeq(:,(44+1):end)=eye(6)*T;
        
        % Form beq matrix of ID loads
        sum_all_c2 = sum(sum_c2,2);
        beq = AllID(t,:)' - sum_all_c2;
        
        % There are no inequality constrains
        A=[];
        b=[];
        
        [x,~,exitflag_aux2] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],opts);
        Activations_ext_aux=[Activations_ext_aux; x'];
        exitflag_aux=[exitflag_aux2; exitflag_aux];
        IDLoadsMatched_aux=[IDLoadsMatched_aux; (Aeq*x + sum_all_c2)'];
    end
end

Activations_ext=[Activations_ext; Activations_ext_aux];
exitflag=[exitflag; exitflag_aux];
IDLoadsMatched=[IDLoadsMatched; IDLoadsMatched_aux];

end


Activations2(1:2:101,:)=Activations_ext(:,1:44);
Activations2(2:2:100,:)=(Activations2(1:2:99,:)+Activations2(3:2:101,:))/2;
Activations=Activations2;

Activations_res=Activations_ext(:,45:end);
Activations_res2(1:2:101,:)=Activations_res;
Activations_res2(2:2:100,:)=(Activations_res2(1:2:99,:)+Activations_res2(3:2:101,:))/2;
Activations_res=Activations_res2;
Reserve_Actuators=Activations_res*T;