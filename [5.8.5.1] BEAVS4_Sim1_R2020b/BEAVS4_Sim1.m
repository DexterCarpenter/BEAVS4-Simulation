%% OSU ESRA 2024
% BEAVS calculations to predict behavior
% Dexter Carpenter
% MATLAB Version R2020b

% In order for this script to run properly, you must have the following
% files in the same folder as the script
    % RocketDataDefault.csv
    % RocketEventDefault.csv
    % RocketVarNames.csv

clear
clc

% setup figures
fig1 = figure(1); figure(fig1); clf

%% Import Data
% The OpenRocket file export must have the following config
    % Must go Edit > Preferences > Units > Default Metric > Ok
    % Select All Variables
    % Format > Field Seperator String > Comma
    % Decimal Places > 3
    % Use Exponential Notation
    % Section Comments > Don't Include Simulation Description
    % Section Comments > Don't Include Flight Events
    % Comment Character > #
    % Keep Field Descriptions

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
if isempty(RocketEvent_file) == true; RocketEvent_file = 'RocketEventDefault.csv'; end
RocketEvent = readtable(RocketEvent_file,'VariableNamingRule','preserve'); clear RocketEvent_file
% Extract values from event data
RocketEvent.Properties.VariableNames = {'Var1','Var2','Name','Var4','Var5','Time','Var7'};
RocketEvent.Time = convertCharsToStrings(RocketEvent.Time);
for i = 1:numel(RocketEvent.Time)
    RocketEvent.Time(i) = regexp(RocketEvent.Time(i),'\d+[\.]?\d*','match');
end
RocketEvent.Time = str2double(RocketEvent.Time);

clc

%% Fill Vars with RocketData

Time = RocketData.Time;
h    = RocketData.Altitude;
V    = RocketData.VerticalVelocity;
Cd   = RocketData.DragCoefficient;

%% Implement Forward Euler with OpenRocket Cd

% extension of BEAVS over TIME (convert to meters)
BEAVSExtn = zeros(numel(RocketData.Time),1) + 6*10^-2;

% Use Forward Euler to calculate velocity and altitude
% Include Cd Scaling
n = 2;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, true);

% Ignore Cd Scaling
n = 3;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, false);

% 110% 
%[Time3, h3, V3, Cd3] = FEuler(RocketData, RocketEvent,  BEAVSExtn*1.1, true);

%% Plot Results

% PLOT
hold on
% plot altitudes
yyaxis left % alt on left
plot(RocketData.Time,h*3.28084)
ylabel('Altitude (ft)');
axis([0 max(max(Time)) 0 12000])
% plot velocity
yyaxis right % velocity on right
plot(RocketData.Time,V*3.28084)
ylabel('Velocity (ft/s)');
axis([0 max(max(Time)) 0 1125.33]);

xlabel('Time (s)');
legend('OpenRocket','Scaled Cd BEAVS','Unscaled Cd BEAVS',...
       'OpenRocket','Scaled Cd BEAVS','Unscaled Cd BEAVS')

%% Summary of Simulation

Name = {'OpenRocket','Scaled Cd BEAVS','Unscaled Cd BEAVS'};
for i = 1:numel(Name)
    Summary(Time(:,i), h(:,i), V(:,i), Cd(:,i), RocketData, Name{:,i})
end

%% Functions

% Determine how much to scale Cd given Fin Extension
% Obtained from OpenRocket Data
% see Cd_Interps.m
function y = ScaleCd(Extsn)
    y = 14.6740.*Extsn + 1;
end

function [Time, h, V, Cd] = FEuler(RocketData, RocketEvent, Extension, ScaleCdBool)
    % Calculate Areas
    BEAVSWdth = 7*10^-2; % width of BEAVS Fin (convert to meters)
    BEAVSFinCnt = 2;     % Number of Fins
    ABeavs = Extension.*BEAVSWdth;
    %ABeavs = ABeavs.*0;
    AFins = 0.00508*0.10795; ATube = 0.02043; ANacelle = (1.25*2.54)*(1/100)^2;
    A = ATube + 4.*AFins + 2.*ANacelle + BEAVSFinCnt.*ABeavs;
    
    % Determine where to start
    BurnoutTime = RocketEvent.Time(4);
    iterStart = find(RocketData.Time==BurnoutTime); % iteration start time
    
    % Allocate vectors for Beavs Data
    Time = RocketData.Time;
    h              = zeros(numel(Time),1);
    h(1:iterStart) = RocketData.Altitude(1:iterStart);
    V              = zeros(numel(Time),1);
    V(1:iterStart) = RocketData.VerticalVelocity(1:iterStart);
    
    % constants
    g = RocketData.GravitationalAcceleration;
    Rd = 287.058;
    rho = (RocketData.AirPressure.*100)./Rd./(RocketData.AirTemperature + 273.15);
    Cd = RocketData.DragCoefficient; Cd(isnan(Cd)) = 0;
    if ScaleCdBool == true
        Cd = ScaleCd(Extension).*Cd;
    end
    m = RocketData.Mass/1000;
    T = RocketData.Thrust;
    
    for i = iterStart:numel(Time)-1
        % Calculate V
        V(i+1) = V(i) + (...
                -g(i)...                                       % gravity
                - 1/2*rho(i)*V(i)*abs(V(i))*Cd(i)*A(i)/m(i)... % drag
                + V(i)/abs(V(i))*T(i)/m(i)...                  % thrust
                )*(Time(i+1) - Time(i));
        % Calculate h
        h(i+1) = h(i) + V(i)*(Time(i+1) - Time(i));
    end
end

function Summary(Time, h, V, Cd, RocketData, Name)
    h = h*3.28084; V = V*3.28084;
    RocketData.Altitude = RocketData.Altitude*3.28084;
    fprintf('\n-----%s-----\n',Name);
    fprintf('Apogee: %0.0f ft\n',max(h));
    fprintf('Apogee Difference: %0.0f ft\n',max(RocketData.Altitude)-max(h));
    fprintf('\n');
%    fprintf('Max Velocity: %0.0f ft/s\n',max(V));
%    fprintf('\n');
    fprintf('Cd_min: %0.2f, Cd_max: %0.2f\n',min(Cd(Cd>0)),max(Cd));
end
