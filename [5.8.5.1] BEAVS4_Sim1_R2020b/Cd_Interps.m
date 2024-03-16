%% Cd Interpolations
% This script is trying to get an understanding of how the Cd changes with
% repsect to the change in BEAVS fin extension
% Dexter Carpenter

clear
clc

fig2 = figure(2); figure(fig2); clf

CdData_file = 'Xx70_Beavs_Drag.csv';
CdData = readtable(CdData_file,'VariableNamingRule','preserve'); clear CdData_file

for k = 1:numel(CdData{1,:})
    for j = 1:numel(CdData{:,1})
        if isnan(CdData{j,k}) == true
            CdData{j,k} = CdData{j-1,k};
        end
    end
end

T  = CdData{:,1};
A = (0:7)*10^-2;
Cd = CdData{:,2:end}';

[T,A] = meshgrid(T,A);
s = surf(A,T,Cd,'EdgeColor','none');
xlabel('BEAVS Fin Extn (m)');
ylabel('Time');
zlabel('Cd');
axis([0 0.07 0 30 0 1.2]);


for i = 2:9
    Cd_prcnt_inc(i-1) = min(CdData{:,i})./min(CdData{:,2});
    fprintf('min: %0.2f, max: %0.2f, %%increase: %0.2f\n',min(CdData{:,i}), max(CdData{:,i}), min(CdData{:,i})./min(CdData{:,2}))
end

Cd_prcnt_inc

coeffs = polyfit(A(:,1), Cd_prcnt_inc, 1)

