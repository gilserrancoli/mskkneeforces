function [AllID, Allmoment_arms]=Reorder_ID_ma(ID,moment_arms,~)

% Author: Gil Serrancolí

%% Reorder experimental data according to dofs

for i=1:length(ID.Coordinates)
    switch ID.Coordinates{i}
        case 'hip_flexion'
            AllID(:,1)=ID.data(:,i+1);
        case 'hip_adduction'
            AllID(:,2)=ID.data(:,i+1);
        case 'hip_rotation'
            AllID(:,3)=ID.data(:,i+1);
        case 'knee_flexion'
            AllID(:,4)=ID.data(:,i+1);
        case 'knee_adduction'
            AllID(:,5)=ID.data(:,i+1);
        case 'knee_ty'
            AllID(:,6)=ID.data(:,i+1);
        case 'ankle_angle'
            AllID(:,7)=ID.data(:,i+1);
        case 'subtalar_angle'
            AllID(:,8)=ID.data(:,i+1);
    end
end

for i=1:length(moment_arms)
    switch moment_arms(i).coordinate
        case 'hip_flexion'
            Allmoment_arms(:,1:44)=moment_arms(i).data(:,2:end);
        case 'hip_adduction'
            Allmoment_arms(:,(44*1+1):44*2)=moment_arms(i).data(:,2:end);
        case 'hip_rotation'
            Allmoment_arms(:,(44*2+1):44*3)=moment_arms(i).data(:,2:end);
        case 'knee_flexion'
            Allmoment_arms(:,(44*3+1):44*4)=moment_arms(i).data(:,2:end);
        case 'knee_adduction'
            Allmoment_arms(:,(44*4+1):44*5)=moment_arms(i).data(:,2:end);
        case 'knee_ty'
            Allmoment_arms(:,(44*5+1):44*6)=moment_arms(i).data(:,2:end);
        case 'ankle_angle'
            Allmoment_arms(:,(44*6+1):44*7)=moment_arms(i).data(:,2:end);
        case 'subtalar_angle'
            Allmoment_arms(:,(44*7+1):44*8)=moment_arms(i).data(:,2:end);
    end
end
