%% StateSpacePlayground
% this script is indented to be run after BEAVS4_Sim5 to use the workspace
% variables. The intention is to mess around with and learn how nonlinear
% state space models work in MATLAB

% This script requires System Identification Toolbox
% https://www.mathworks.com/help/ident/index.html?searchHighlight=System%20Identifying%20Toolbox&s_tid=srchtitle_support_results_1_System%20Identifying%20Toolbox

% Measured input data set
U = Cd(:,4);
% Measured output data set
Y = V(:,4);
% Create Neural State-Space Object
% Singe Input, Single Output (SISO)

tt = [U Y];
nx = 1;
Ts = 0.01;
ct_sys = ssest(tt, nx,'Ts',Ts);

pidTuner(ct_sys, 'PID');














