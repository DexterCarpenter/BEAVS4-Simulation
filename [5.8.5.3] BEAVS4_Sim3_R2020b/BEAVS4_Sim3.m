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

% Fill Vars with RocketData to start (these values will be used as
% reference
Time = RocketData.Time;
h    = RocketData.Altitude;
V    = RocketData.VerticalVelocity;
Cd   = RocketData.DragCoefficient;

%% Additional Script Inputs

BladeWdth = 7*10^-2; % width of BEAVS Blade (convert to meters)
BladeCnt = 2;        % Number of Blades
BladeExtnRate = 0.074;   % meters per second

% extension of BEAVS over TIME (convert to meters)
BEAVSExtn = zeros(numel(RocketData.Time),1) + 6*10^-2;
BEAVSExtnMAX = 0.06; % max extension of BEAVS in m
BurnoutTime = RocketEvent.Time(4);
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time
BEAVSExtn(1:iterStart) = 0;
for i = iterStart:numel(BEAVSExtn)
    if BEAVSExtn(i) >= BEAVSExtnMAX
        break
    end
    BEAVSExtn(i+1) = BEAVSExtn(i) + BladeExtnRate*( Time(i+1) - Time(i) );
end

%% Implement Forward Euler with OpenRocket Cd

% Use Forward Euler to calculate velocity and altitude

% Use OpenRocket Cd
n = 2;
Cd_rocket = RocketData.DragCoefficient;
Cd_BEAVS = 1.10;
FeedbackBool = false;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt, FeedbackBool);

% Low End Cd Estimation (Least Beneficial for BEAVS)
n = 3;
Cd_rocket = min(RocketData.DragCoefficient);
Cd_BEAVS = 1.00;
FeedbackBool = false;
[Time(:,n), h(:,n), V(:,n), Cd(:,n)] = FEuler(RocketData, RocketEvent, BEAVSExtn, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt, FeedbackBool);

%% Plot Results

% PLOT
figure(fig1);
hold on
% plot altitudes
yyaxis left % alt on left
plot(RocketData.Time,h*3.28084)
ylabel('Altitude (ft)');
lims = [0 max(max(Time)) 0 12000];
lims(2) = 50;
axis(lims);
idx = find(RocketData.Time==RocketEvent.Time(4)); % burnout start time idx
fill([Time(idx:end,1)' flip(Time(idx:end,2)')], [h(idx:end,1)' flip(h(idx:end,2)')]*3.28084,...
    [0 0.4470 0.7410], 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% plot velocity
yyaxis right % velocity on right
plot(RocketData.Time,V*3.28084)
ylabel('Velocity (ft/s)');
lims = [0 max(max(Time)) 0 1125.33];
lims(2) = 50;
axis(lims);
xlabel('Time (s)');
xline(RocketEvent.Time(4),'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
yyaxis left
yline(10000,'--','Target Apogee','LabelVerticalAlignment','middle','LabelHorizontalAlignment','right');
legend('OpenRocket','Max Cd BEAVS','Min Cd BEAVS','Location','southeast')
title('BEAVS 4.0 Control Range Graphic');
grid on


figure(fig2);
hold on
plot(Time(:,1:2),Cd(:,1:2));
lims = [0 50 0 1.2];
axis(lims);
xline(RocketEvent.Time(4),'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
xlabel('Time (s)');
ylabel('Cd total');
title('Cd over Time');
legend('w/o BEAVS','w/ BEAVS','Location','SouthEast');
grid on

%% Summary of Simulation

Name = {'OpenRocket','Scaled Cd BEAVS','Unscaled Cd BEAVS'};
for i = 1:n
    Summary(Time(:,i), h(:,i), V(:,i), Cd(:,i), RocketData, Name{:,i})
end

%% Functions

% Determine how much to scale Cd given Fin Extension
% Obtained from OpenRocket Data
% see Cd_Interps.m
function Cd_total = ScaleCd(Cd_rocket, Cd_i, Ai, Aref)
    %y = 14.6740.*Extsn + 1;
    Cd_total = Cd_rocket + Cd_i*(Ai/Aref);
end

function [Time, h, V, Cd] = FEuler(RocketData, RocketEvent, Extension, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt, FeedbackBool)
  % This function takes OpenRocket Data, BEAVS blade extension, and Cd's
  % and applies Forward Euler to calculate altitude and velocity over time
% INPUTS ------------------------------------------------------------------
    % RocketData   | Table from OpenRocket
    % RocketEvent  | Table of events from OpenRocket
    % Extension    | vector containing BEAVS blade extension over the
    %                OpenRocket Time vector
    % Cd_rocket    | desired coefficient of drag for the rocket w/o BEAVS
    %                can be vector w/ respect to OpenRocket Time
    % Cd_BEAVS     | desired coefficient of drag for BEAVS baldes ALONE
    %                can be vector w/ respect to OpenRocket Time
    % FeedbackBool | true - use getfeedback.m
    %              | false - do not use getfeedback.m
% OUTPUTS -----------------------------------------------------------------
    % Time | time vector, the same as OpenRocket time vector, just helpful
    % h    | altitude vector w/ respect to Time
    % V    | velocity vector w/ respect to Time
    % Cd   | calculated Cd vector for reference
    
    % if feedback is not used, pre-calculate area vector
    if FeedbackBool == false
        % Reference Area
        A.RefOpen = RocketData.ReferenceArea.*(0.01).^2; % convert to m

        % Calculate BEAVS Area
        A.BEAVS = Extension.*BladeWdth.*BladeCnt; % calculate BEAVS area

        % Unused
        %A.Fins = 0.00508*0.10795;
        %A.Tube = 0.02043;
        %A.Nacelle = (1.25*2.54)*(1/100)^2;
        %A.RefBEAVS = A.Tube + A.BEAVS;

        % Calculate Cd TOTAL
        %Cd = RocketData.DragCoefficient; Cd(isnan(Cd)) = 0;
        Cd = Cd_rocket + Cd_BEAVS.*(A.BEAVS./A.RefOpen);
    end
    
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
        if FeedbackBool == true
            getfeedback([Time(i) Time(i+1)]);
        end
        
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
