%% OSU ESRA 2024
% BEAVS calculations to predict behavior
% Dexter Carpenter

% The purpose of this script is to provide an idea of what a theoretical
% maximum air braking power BEAVS 4.0 could achieve.
% Since many changes are being made to the rocket as time progresses,
% having a script to quickly see how apogee changes with different motors,
% Cds, BEAVS blade sizes, etc. it quickly becomes too much to handle
% otherwise

% This script takes OpenRocket Data exports along with a set of additional
% BEAVS variables to calculate apogees

% In order for this script to run properly, you must have the following
% files in the same folder as the script
    % RocketDataDefault.csv
    % RocketEventDefault.csv
    % RocketVarNames.csv

clear
clc

% setup figures
fig1 = figure(1); figure(fig1); clf
fig2 = figure(2); figure(fig2); clf
fig3 = figure(3); figure(fig3); clf

%% Import Data
% OpenRocket Data Export config:
    % Must go Edit > Preferences > Units > Default Metric > Ok
    % Select All Variables
    % Format > Field Seperator String > Comma
    % Decimal Places > 3
    % Use Exponential Notation
    % Section Comments > Don't Include Simulation Description
    % Section Comments > Don't Include Flight Events
    % Comment Character > #
    % Keep Field Descriptions
    
% OpenRocket Event Export config:
    % Variables to Export > None
    % Format Settings > csv
    % Include Flight Events
    % Comment Character > #

% Get OpenRocket Data
disp('You can Press ENTER to select default data');
RocketData_file = input('Rocket Data File Name: ',"s");
if isempty(RocketData_file) == true; RocketData_file = 'RocketDataDefault.csv'; fprintf('Default Selected\n'); end
RocketData = readtable(RocketData_file,'VariableNamingRule','preserve'); clear RocketData_file

% Set Table Variable Names from csv
RocketVarNames = readtable('RocketVarNames.csv','VariableNamingRule','preserve');
RocketData.Properties.VariableNames = RocketVarNames.Properties.VariableNames;

% Get OpenRocket Event Data
RocketEvent_file = input('Rocket Event File Name: ',"s");
if isempty(RocketEvent_file) == true; RocketEvent_file = 'RocketEventDefault.csv'; fprintf('Default Selected\n'); end
RocketEvent = readtable(RocketEvent_file,'VariableNamingRule','preserve'); clear RocketEvent_file
% Extract values from event data
RocketEvent.Properties.VariableNames = {'Var1','Var2','Name','Var4','Var5','Time','Var7'};
RocketEvent.Time = convertCharsToStrings(RocketEvent.Time);
for i = 1:numel(RocketEvent.Time)
    RocketEvent.Time(i) = regexp(RocketEvent.Time(i),'\d+[\.]?\d*','match');
end
RocketEvent.Time = str2double(RocketEvent.Time);

% fill NaN values
RocketData = fillmissing(RocketData, 'previous');

%% Additional Script Inputs

BladeWdth = 7*10^-2; % width of BEAVS Blade (convert to meters)
BladeCnt = 2;        % Number of Blades
BladeExtnRate = 110e-3;   % meters per second
BEAVSExtnMAX = 0.06; % max extension of BEAVS in m

% extension of BEAVS over TIME (convert to meters)
BEAVSExtn = zeros(numel(RocketData.Time),1) + 6*10^-2;
BurnoutTime = RocketEvent.Time(4); % define time of burnout
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time
BEAVSExtn(1:iterStart) = 0; % no extension in boost phase
% ramp up BEAVS extension until fully extended at extension rate
for i = iterStart:numel(BEAVSExtn)
    if BEAVSExtn(i) >= BEAVSExtnMAX; break; end
    BEAVSExtn(i+1) = BEAVSExtn(i) + BladeExtnRate*( RocketData.Time(i+1) - RocketData.Time(i) );
end

%% Implement Forward Euler with OpenRocket Cd
% Use Forward Euler to calculate velocity and altitude

% OpenRocket Simulation, No BEAVS
n = 1;
SimName(n) = {'OpenRocket'};
Time = RocketData.Time;
h    = RocketData.Altitude;
V    = RocketData.VerticalVelocity;
Cd   = RocketData.DragCoefficient;

% Include BEAVS, Constant Cd
% n = 2;
% SimName(n) = {'Constant Cd'};
% Cd_rocket = RocketData.DragCoefficient;
% Cd_BEAVS = 1.10;
% [Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt);

% Include BEAVS, Vary Cd
n = 2;
SimName(n) = {'Varying Cd'};
Cd_rocket = RocketData.DragCoefficient;
Aref = RocketData.ReferenceArea.*(0.01).^2; % convert to m
Cd_BEAVS = InterpCd(BEAVSExtn.*BladeWdth,Aref);
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt);

% Include BEAVS, Feedback & Constant Cd
n = 3;
SimName(n) = {'Feedback'};
Cd_rocket = RocketData.DragCoefficient;
[Time(:,n), h(:,n), V(:,n), Cd(:,n), Extn, ExtnDesire] = FEulerFeedback(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BEAVSExtnMAX);

% Runge-Kutta (ode45)
% n = 4;
% SimName(n) = {'Runge-Kutta'};
% Cd_rocket = RocketData.DragCoefficient;
% Cd_BEAVS = 1.10;
% FeedbackBool = false;
% [Time(:,n), h(:,n), V(:,n), Cd(:,n)] = RKutta(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt);

%% Figure 1
% Altitude and Cd

figure(fig1);
hold on
title('BEAVS 4.0 Control Range Graphic');
xlabel('Time (s)');

% LEFT --------------------------------------------------------------------
yyaxis left % Altitude
plot(Time(:,1:n),h*3.28084);
ylabel('Altitude (ft)');
lims = [0 50 0 12000]; axis(lims);

% Target Apogee Line
yline(10000,'--','Target Apogee','LabelVerticalAlignment','middle','LabelHorizontalAlignment','right');

% Shade control region
idx = find(RocketData.Time==BurnoutTime); % burnout start time idx
fill([Time(idx:end,1)' flip(Time(idx:end,2)')], [h(idx:end,1)' flip(h(idx:end,2)')]*3.28084,...
    [0 0.4470 0.7410], 'FaceAlpha', 0.05, 'EdgeColor', 'none');

% RIGHT -------------------------------------------------------------------
yyaxis right % Coefficient of Drag
plot(Time(:,1:n),Cd(:,1:n));
ylabel('Cd total');
lims = [0 50 0 1.8]; axis(lims);

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','southeast')
grid on

%% Figure 2
% Velocity

figure(fig2);
hold on
title('Cd over Time');
xlabel('Time (s)');

% Velocity
plot(Time(:,1:n),V*3.28084);
ylabel('Velocity (ft/s)');
lims = [0 50 0 1125.33]; axis(lims);

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','SouthEast');
grid on

%% Figure 3

figure(fig3);
hold on
title('Extension & Desired Extension');
xlabel('Time (s)');

plot(Time(:,3),[Extn ExtnDesire]*1000);
lims = [0 25 0 BEAVSExtnMAX*1000]*1.1; axis(lims);
ylabel('Extension (mm)');

yline(BEAVSExtnMAX*1000,'--','Max Extension');
legend('Actual Extn','Desired Extn','Location','SouthEast');
grid on

%% Summary of Simulation

% output summary to cmd window for each simulation
for i = 1:n
    Summary(Time(:,i), h(:,i), V(:,i), Cd(:,i), RocketData, RocketEvent, SimName{:,i})
end
