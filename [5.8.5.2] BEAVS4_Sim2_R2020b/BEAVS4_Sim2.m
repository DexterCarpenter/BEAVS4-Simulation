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

% Include BEAVS, Vary Cd
n = 2;
SimName(n) = {'Varying Cd'};
Cd_rocket = RocketData.DragCoefficient;
Cd_BEAVS = 1.10;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt);

% Include BEAVS, Minimum Cd
n = 3;
SimName(n) = {'Minimum Cd'};
Cd_rocket = min(RocketData.DragCoefficient);
Cd_BEAVS = 1.10;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt);

%% Figure 1
% Apogee and Cd

figure(fig1);
hold on
title('BEAVS 4.0 Control Range Graphic');
xlabel('Time (s)');

% LEFT --------------------------------------------------------------------
yyaxis left % Altitude
plot(RocketData.Time,h*3.28084);
ylabel('Altitude (ft)');
lims = [0 50 0 12000]; axis(lims);

% Target Apogee Line
yline(10000,'--','Target Apogee','LabelVerticalAlignment','middle','LabelHorizontalAlignment','right');

% Shade control region
idx = find(RocketData.Time==RocketEvent.Time(4)); % burnout start time idx
fill([Time(idx:end,1)' flip(Time(idx:end,2)')], [h(idx:end,1)' flip(h(idx:end,2)')]*3.28084,...
    [0 0.4470 0.7410], 'FaceAlpha', 0.05, 'EdgeColor', 'none');

% RIGHT -------------------------------------------------------------------
yyaxis right % Coefficient of Drag
plot(Time(:,1:3),Cd(:,1:3));
ylabel('Cd total');
lims = [0 50 0 1.5]; axis(lims);

% Motor Cutoff and Legend
xline(RocketEvent.Time(4),'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','southeast')
grid on

%% Figure 2
% Velocity

figure(fig2);
hold on
title('Cd over Time');
xlabel('Time (s)');

% Velocity
plot(RocketData.Time,V*3.28084);
ylabel('Velocity (ft/s)');
lims = [0 50 0 1125.33]; axis(lims);

% Motor Cutoff and Legend
xline(RocketEvent.Time(4),'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','SouthEast');
grid on

%% Summary of Simulation

% output summary to cmd window for each simulation
for i = 1:n
    Summary(Time(:,i), h(:,i), V(:,i), Cd(:,i), RocketData, RocketEvent, SimName{:,i})
end

%% Functions --------------------------------------------------------------

function [Time, h, V, Cd] = FEuler(RocketData, RocketEvent, Extension, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt)
  % This function takes OpenRocket Data, BEAVS blade extension, and Cd's
  % and applies Forward Euler to calculate altitude and velocity over time
% INPUTS ------------------------------------------------------------------
    % RocketData  | Table from OpenRocket
    % RocketEvent | Table of events from OpenRocket
    % Extension   | vector containing BEAVS blade extension over the
    %               OpenRocket Time vector
    % Cd_rocket   | desired coefficient of drag for the rocket w/o BEAVS
    %               can be vector w/ respect to OpenRocket Time
    % Cd_BEAVS    | desired coefficient of drag for BEAVS baldes ALONE
    %               can be vector w/ respect to OpenRocket Time
    % BladeWdth   | Width of BEAVS blades
    % BladeCnt    | Number of BEAVS baldes
% OUTPUTS -----------------------------------------------------------------
    % Time | time vector, the same as OpenRocket time vector, just helpful
    % h    | altitude vector w/ respect to Time
    % V    | velocity vector w/ respect to Time
    % Cd   | calculated Cd vector for reference
    
    % Reference Area
    A.RefOpen = RocketData.ReferenceArea.*(0.01).^2; % convert to m
    
    % Calculate BEAVS Area
    A.BEAVS = Extension.*BladeWdth.*BladeCnt; % calculate BEAVS area
    
    % Calculate Cd TOTAL
    Cd = Cd_rocket + Cd_BEAVS.*(A.BEAVS./A.RefOpen);
    
    % Determine where to begin iteration (start at burnout)
    BurnoutTime = RocketEvent.Time(4);
    iterStart = find(RocketData.Time==BurnoutTime); % iteration start time
    
    % Allocate vectors for Beavs Data
    Time           = RocketData.Time;
    h              = zeros(numel(Time),1);             % total vector
    h(1:iterStart) = RocketData.Altitude(1:iterStart); % boost phase
    V              = zeros(numel(Time),1);                     % total vector
    V(1:iterStart) = RocketData.VerticalVelocity(1:iterStart); % boost phase
    
    % constants
    g = RocketData.GravitationalAcceleration;
    Rd = 287.058;
    rho = (RocketData.AirPressure.*100)./Rd./(RocketData.AirTemperature + 273.15);
    m = RocketData.Mass/1000;
    T = RocketData.Thrust;
    
    for i = iterStart:numel(Time)-1
        % Calculate V
        V(i+1) = V(i) + (...
                -g(i)...                                               % gravity
                - 1/2*rho(i)*V(i)*abs(V(i))*Cd(i)*A.RefOpen(i)/m(i)... % drag
                + V(i)/abs(V(i))*T(i)/m(i)...                          % thrust
                )*(Time(i+1) - Time(i));
        % Calculate h
        h(i+1) = h(i) + V(i)*(Time(i+1) - Time(i));
    end
end

% Format and Summarize Data to Command Window
function Summary(Time, h, V, Cd, RocketData, RocketEvent, Name)
    % Convert to ft
    h = h*3.28084; V = V*3.28084;
    RocketData.Altitude = RocketData.Altitude*3.28084;
    % find coast phase start time
    iterStart = find(RocketData.Time==RocketEvent.Time(4));
    % print to cmd window
    fprintf('\n----- %s -----\n',Name);
    fprintf('Apogee: %0.0f ft\n',max(h));
    fprintf('Apogee Difference: %0.0f ft\n',max(RocketData.Altitude)-max(h));
    fprintf('\n');
    fprintf('Max Velocity: %0.0f ft/s\n',max(V));
    fprintf('\n');
    fprintf('Cd_min: %0.2f, Cd_max: %0.2f\n',min(Cd(iterStart:end)),max(Cd(iterStart:end)));
end
