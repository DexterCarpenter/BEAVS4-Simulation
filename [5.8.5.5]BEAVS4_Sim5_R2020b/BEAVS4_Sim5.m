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
fig4 = figure(4); figure(fig4); clf
fig5 = figure(5); figure(fig5); clf
fig6 = figure(6); figure(fig6); clf
fig7 = figure(7); figure(fig7); clf

%% Import Simulation Data -------------------------------------------------
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
fprintf('You can Press ENTER to select default data\n');
fprintf('You can enter "1" to select data set 1\n');
fprintf('Example: SomeWildRocketData.csv\n');
RocketData_file = input('Rocket Data File Name: ',"s");
if (isempty(RocketData_file) == true) || RocketData_file == "0" % SELECT DEFAULT DATA SET
    RocketData_file = 'RocketDataSet0.csv'; fprintf('Default Selected\n');
    RocketEvent_file = 'RocketEventSet0.csv'; fprintf('Rocket Event File Name: \nDefault Selected\n');
elseif RocketData_file == "1" % SELECT FIRST DATA SET
    RocketData_file = 'RocketDataSet1.csv'; fprintf('Data Set 1 Selected\n');
    RocketEvent_file = 'RocketEventSet1.csv'; fprintf('Rocket Event File Name: \nData Set 1 Selected\n');
elseif RocketData_file == "2" % SELECT APRIL BROTHERS DATA SET
    RocketData_file = 'RocketDataSet2.csv'; fprintf('Data Set 2 Selected (April Brothers Set)\n');
    RocketEvent_file = 'RocketEventSet2.csv'; fprintf('Rocket Event File Name: \nData Set 2 Selected (April Brothers Set)\n');
else % CUSTOM DATA SET
    RocketEvent_file = input('Rocket Event File Name: ',"s");
end

% Extract Rocket DATA
RocketData = readtable(RocketData_file,'VariableNamingRule','preserve'); clear RocketData_file

% Set Table Variable Names from csv
RocketVarNames = readtable('RocketVarNames.csv','VariableNamingRule','preserve');
RocketData.Properties.VariableNames = RocketVarNames.Properties.VariableNames;

% Extract Rocket EVENTS
RocketEvent = readtable(RocketEvent_file,'VariableNamingRule','preserve'); clear RocketEvent_file
RocketEvent.Properties.VariableNames = {'Var1','Var2','Name','Var4','Var5','Time','Var7'};
RocketEvent.Time = convertCharsToStrings(RocketEvent.Time);
for i = 1:numel(RocketEvent.Time)
    RocketEvent.Time(i) = regexp(RocketEvent.Time(i),'\d+[\.]?\d*','match');
end
RocketEvent.Time = str2double(RocketEvent.Time);

% fill NaN values with previous value
RocketData = fillmissing(RocketData, 'previous');

%% Import Actual Data -----------------------------------------------------

% There were two flight computers in flight for the launch: 11445 and 11439
% which will be reffered to as "1" and "2" respectively. 

% Import data from April Brothers 2024 Launch
AprilBrothers1 = readtable("April_Brothers_TeleMetrum_11445.csv",'VariableNamingRule','preserve'); % TeleMetrum 1
AprilBrothers2 = readtable("April_Brothers_TeleMetrum_11439.csv",'VariableNamingRule','preserve'); % TeleMetrum 2

% Find Motor Cutoff
AprilBrothersCoastTable1 = AprilBrothers1(:,[5 7 11 12]); % extract height, time, speed
coastRows = strcmp(AprilBrothers1.state_name, 'coast') == 0; % logical array for not 'coast'
AprilBrothersCoastTable1(coastRows,:) = []; % delete non-coast rows
AprilBrothersCutoff1 = AprilBrothersCoastTable1.time(1); % extract time of cutoff

AprilBrothersCoastTable2 = AprilBrothers2(:,[5 9 13 14]); % extract height, time, speed
coastRows = strcmp(AprilBrothers2.state_name, 'coast') == 0; % logical array for not 'coast'
AprilBrothersCoastTable2(coastRows,:) = []; % delete non-coast rows
AprilBrothersCutoff2 = AprilBrothersCoastTable2.time(1); % extract time of cutoff

% Remove repeated values
[~,ia] = unique(AprilBrothersCoastTable1.height);
LookupTable1 = AprilBrothersCoastTable1(ia,:);
[~,ia] = unique(AprilBrothersCoastTable2.height);
LookupTable2 = AprilBrothersCoastTable2(ia,:);

% Create Velocity Lookup Table
% velocity lookup table gives the PID a target velocity to achieve based
% upon a given altitude. the idea is to have the PID system control the
% velocity instead of a height. Version 4 used a target altitude of
% 10,000ft instead of a velocity lookup table
LookupRes = 0.01;
VelLookup(:,2) = (0:LookupRes:max(LookupTable1.speed*1.05))';
VelLookup(:,1) = polyval(    polyfit(LookupTable1.speed,LookupTable1.height,5) , VelLookup(:,2)    );

figure(fig7); hold on; grid on;
plot(VelLookup(:,2),VelLookup(:,1));
plot(LookupTable1.speed,LookupTable1.height);
plot(LookupTable2.speed,LookupTable2.height);
title('Altitude vs. Velocity Lookup Table');
xlabel('Velocity (ft/s)'); ylabel('Altitude (ft)');

%% Additional Script Inputs -----------------------------------------------

% Determine coast phase begin time and index number
BurnoutTime = RocketEvent.Time(4); % define time of burnout
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time

% BEAVS Characteristic Parameters
BladeWdth = 7*10^-2;    % width of BEAVS Blade (convert to meters)
BladeCnt = 2;           % Number of Blades
BladeExtnRate = 110e-3; % meters per second
BladeExtnMAX = 0.06;    % max extension of BEAVS in m

% Create BladeExtn Vector
BladeExtn = zeros(numel(RocketData.Time),1) + 6*10^-2; % extension of BEAVS over TIME (convert to meters)
BladeExtn(1:iterStart) = 0; % no extension in boost phase
% ramp up BEAVS extension until fully extended at extension rate
for i = iterStart:numel(BladeExtn)-1
    % if max extension is hit, stop iterating
    if BladeExtn(i) >= BladeExtnMAX; break; end
    % Extend blades
    BladeExtn(i+1) = BladeExtn(i) + BladeExtnRate*( RocketData.Time(i+1) - RocketData.Time(i) );
end

% NOISE
% degree of noise to use in simulations
% 0-100 corrolating to a percent
% see Noise.m
NoiseDeg = 5;

%% SIMULATIONS ------------------------------------------------------------
% Use Forward Euler to calculate velocity and altitude

% various constants
Cd_rocket = RocketData.DragCoefficient;
Cd   = RocketData.DragCoefficient;
hTarg = 10000*0.3048;
A_ref = RocketData.ReferenceArea(1).*(0.01).^2; % convert to m2
rho = (RocketData.AirPressure.*100)./287.058./(RocketData.AirTemperature + 273.15); % convert to kg/m3rho = 

% Pure OpenRocket Simulation, No BEAVS
n = 1;
SimName(n) = {'OpenRocket'};
Time = RocketData.Time;
h    = RocketData.Altitude;
V    = RocketData.VerticalVelocity;

% Maximum Braking
n = n+1;
SimName(n) = {'Maximum Braking'};
[Time(:,n), h(:,n), V(:,n), Cd(:,n), Fb(:,n), Fn(:,n)] = ...
    FEuler(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX,NoiseDeg,"No Feedback",BladeExtn);
Fb(700:end,:) = 0;
Fn(700:end,:) = 0;

% Prediction
n = n+1;
SimName(n) = {'Prediction'};
[Time(:,n), h(:,n), V(:,n), Cd(:,n), PredBladeExtn, PredBladeExtnDesire] = ...
    FEuler(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX,NoiseDeg,"Pred");

% PID
n = n+1;
SimName(n) = {'PID'};
Kp = 1.200e-04; % Kp = 1.320e-04;
Ki = 1.000e-09; % Ki = 1.000e-07;
Kd = 6.500e-05; % Kd = 1.004e-06;
[Time(:,n), h(:,n), V(:,n), Cd(:,n), PidBladeExtn, PIDu, err, Vtarg] = ...
    FEuler(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX,NoiseDeg,"PID",VelLookup,Kp,Ki,Kd);

%% Figure 1
% Altitude and Cd

figure(fig1);
hold on
title('BEAVS 4.0 Control Range Graphic');
xlabel('Time (s)');

% LEFT AXIS ---------------------------------------------------------------
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

% RIGHT AXIS --------------------------------------------------------------
yyaxis right % Coefficient of Drag
plot(Time(:,1:n),Cd(:,1:n));
ylabel('Cd total');
lims = [0 50 0 2]; axis(lims);

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','southeast')
grid on

%% Figure 2
% Velocity

figure(fig2);
hold on
title('Velocity vs. Time');
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

plot(Time(:,1),[PredBladeExtn PredBladeExtnDesire]*1000);
lims = [0 25 0 BladeExtnMAX*1000]*1.1; axis(lims);
ylabel('Extension (mm)');

yline(BladeExtnMAX*1000,'--','Max Extension');
legend('Actual Extn','Desired Extn','Location','SouthEast');
grid on

%% Figure 4

figure(fig4)
hold on
title('PID Graphic');
xlabel('Time (s)');

plot(Time(:,1),[Vtarg V(:,4) err]*3.28084);
lims = [0 30 -200 1000]; axis(lims);
ylabel('Velocity (ft/s)');

TimeTrgtRched = Time(h(:,4)>=3048,4);
if isempty(TimeTrgtRched)
    TimeTrgtRched = Time(h(:,4)==max(h(:,4)),4);
else
    TimeTrgtRched = TimeTrgtRched(1);
end
xline(TimeTrgtRched,'--',num2str(max(h(:,4))*3.281),'LabelHorizontalAlignment','left');

legend('Target','Actual','Error','Location','NorthEast');
grid on

%% Figure 5

figure(fig5);
hold on

title('Forces During Coast Phase');
xlabel('Time (s)');

plot(Time(:,1),abs(Fb(:,2))); % Force on singular BEAVS blade
plot(Time(:,1),abs(Fn(:,2))); % Force on coupler
lims = [0 40 0 2000]; axis(lims);
ylabel('Force (N)');

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend({'Blade Force', 'Coupler Force'},'Location','NorthEast');
grid on

BladeExtnALLOW = InterpCdInverse(2.*abs(Fb(:,2))./rho./A_ref./(V(:,2).^2),A_ref)./BladeWdth.*1000; % convert to mm

%% Figure 6

% Plot
figure(fig6);
hold on

title('Compare to Real');
xlabel('Time (s)');

% LEFT AXIS -----------------------------------------------------------
yyaxis left % Altitude
plot(Time, h*3.281);
plot(AprilBrothers1.time, AprilBrothers1.height*3.281,'-','LineWidth',2);
plot(AprilBrothers2.time, AprilBrothers2.height*3.281,'--','LineWidth',2);

ylabel('Altitude (ft)');
lims = [0 50 0 12000]; axis(lims);

% Target Apogee Line
yline(10000,'--','Target Apogee','LabelVerticalAlignment','middle','LabelHorizontalAlignment','right');

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Sim Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','right');
xline(mean([AprilBrothersCutoff1 AprilBrothersCutoff2]),'-',{'Real Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');

grid on

% LEFT AXIS -----------------------------------------------------------
yyaxis right % Velocity
plot(Time, V*3.281);
plot(AprilBrothers1.time, AprilBrothers1.speed*3.281,'-','LineWidth',2);
plot(AprilBrothers2.time, AprilBrothers2.speed*3.281,'--','LineWidth',2);

ylabel('Velocity (ft/s)');
lims = [0 50 -1000 1500]; axis(lims);

legend([SimName, 'April Brothers 1', 'April Brothers 2'],'Location','SouthEast');

%% Summary of Simulation

% output summary to cmd window for each simulation
for i = 1:n
    Summary(Time(:,i), h(:,i), V(:,i), Cd(:,i), RocketData, RocketEvent, SimName{:,i})
end

figure(fig5);