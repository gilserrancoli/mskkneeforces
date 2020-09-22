function stop=IterationPlots(x,optimvalues,state,params_mod,params,lMT1,lMt2,lMt3,vMT1,vMT2,vMT3,nmuscles,Allmoment_arms1,Allmoment_arms2,Allmoment_arms3,AllID1,AllID2,AllID3,pos_nonzero_ma,NC_initial_og1,NC_initial_og2,NC_initial_og3,SV_initial,pos_musclesWithSynergies,pos_musclesWithNoSynergies,muscles,muscleswithEMG_reliable,EMG_og1,EMG_og2,EMG_og3,Contact_info,ID_noCF1,ID_noCF2,ID_noCF3,Options)

% Author: Gil Serrancolí

%% Save values of the design variables at each iteration

if ~exist('iterResults.mat')
    xRes(1).x=x;
    save('iterResults','xRes')
else
    vars = whos('-file','iterResults.mat');
    load iterResults
    l=length(xRes);
    xRes(l+1).x=x;
    save('iterResults','xRes');
end
stop=false;



