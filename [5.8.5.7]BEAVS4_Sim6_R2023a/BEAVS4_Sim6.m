%% OSU ESRA 2024
% BEAVS calculations to predict behavior
% Dexter Carpenter

% VERSION 6
% Version 6 is intended to be use by outside sources and is updated for
% MATLAB R2023a

% The purpose of this script is to provide an idea of what a theoretical
% maximum air braking power BEAVS 4.0 could achieve.
% Since many changes are being made to the rocket as time progresses,
% having a script to quickly see how apogee changes with different motors,
% Cds, BEAVS blade sizes, etc. it quickly becomes too much to handle
% otherwise

% This script takes OpenRocket Data exports along with a set of additional
% BEAVS variables to calculate apogees

%% How to Import Simulation Data

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

%% -------------------------------- BEGIN ---------------------------------

clear
clc

%% Check Version & Add-ons

if version('-release') ~= "2023a"; warning('Please use MATLAB R2023a'); end

try % make sure 'Statistics and Machine Learning Toolbox' is installed
    matlab.addons.isAddonEnabled('Statistics and Machine Learning Toolbox');
catch
    warning('Please install "Statistics and Machine Learning Toolbox" Add-On');
end

%% Add folders to PATH

% add current folder and subfolders to path
addpath(genpath(pwd));

%% Setup Figures

% make is such that figures are automatically docked
set(0,'DefaultFigureWindowStyle','docked')

% explicitly define figures and clear them
fig1 = figure(1); figure(fig1); clf
fig2 = figure(2); figure(fig2); clf
fig3 = figure(3); figure(fig3); clf
fig4 = figure(4); figure(fig4); clf
fig5 = figure(5); figure(fig5); clf
fig6 = figure(6); figure(fig6); clf
fig7 = figure(7); figure(fig7); clf
    
%% Import Data

fprintf('------------------------- BEAVS 4.0 SIMULATION 6 --------------------------\n');

% Get OpenRocket Data
fprintf('\nChoose Data Source\n');
fprintf('    Press "ENTER" to select default data\n');
fprintf('    Input "1" to select SPACEPORT AMERICA 2024 data set\n');
fprintf('    Input "2" to select APRIL 2024 BROTHERS data set\n');
fprintf('    Else, input a specific file name to choose that file\n');
fprintf('    Example: SomeWildRocketData.csv\n');
RocketData_file = input('->  Rocket Data File Name: ',"s");

switch RocketData_file
    case {"","0"} % SELECT DEFAULT DATA SET
        RocketData_file  = 'DEFAULT_DataSet.csv';  fprintf('Default Data Selected\n');
        RocketEvent_file = 'DEFAULT_EventSet.csv'; fprintf('Default Event Selected\n');
    case "1" % SELECT SPACEPORT 2024 SIMULATION
        RocketData_file  = 'SACup_DataSet.csv';  fprintf('Data Set 1 Selected (Spaceport 2024 Set)\n');
        RocketEvent_file = 'SACup_EventSet.csv'; fprintf('Event Set 1 Selected (Spaceport 2024 Set)\n');
    case "2" % SELECT APRIL BROTHERS DATA SET
        RocketData_file  = 'AprilBrothers_DataSet.csv';  fprintf('Data Set 2 Selected (April Brothers Set)\n');
        RocketEvent_file = 'AprilBrothers_EventSet.csv'; fprintf('Event Set 2 Selected (April Brothers Set)\n');
    otherwise % CUSTOM DATA SET
        RocketEvent_file = input('->  Rocket Event File Name: ',"s");
end

fprintf('\n---------------------------------------------------------------------------\n');
fprintf('Extracting Rocket Data...\n');
tic

% Extract Rocket DATA
RocketData = readtable(RocketData_file,'VariableNamingRule','preserve'); clear RocketData_file;

% Set Table Variable Names from csv
RocketVarNames = readtable('RocketVarNames.csv','VariableNamingRule','preserve');
RocketData.Properties.VariableNames = RocketVarNames.Properties.VariableNames;

% Extract Rocket EVENTS
RocketEvent = readtable(RocketEvent_file,'VariableNamingRule','preserve'); clear RocketEvent_file;
RocketEvent.Properties.VariableNames = {'Var1','Var2','Name','Var4','Var5','Time','Var7'};
RocketEvent.Time = convertCharsToStrings(RocketEvent.Time);

% reformat data so it's readable
for i = 1:numel(RocketEvent.Time)
    % when imported from OpenRocket, time data is of the format
        % "t=0.00" or "t=263.561"
    % this deletes the "t=" portion leaving only the numerical part
    RocketEvent.Time(i) = regexp(RocketEvent.Time(i),'\d+[\.]?\d*','match');
end
% convert the data type from strings to doubles
RocketEvent.Time = str2double(RocketEvent.Time);

% fill NaN values with previous value to avoid errors
RocketData = fillmissing(RocketData, 'previous');

fprintf('\nExtracting Rocket Data Complete!\n');
toc

%% Import Velocity Lookup Table Reference Data ----------------------------

% There were two flight computers in flight for the launch: 11445 and 11439
% which will be reffered to as "1" and "2" respectively.

% Flight computer 1 seems much more reliable than 2, therefore it is used.

% Ask the user what velocity curve they wish to use for the simulation

fprintf('\n---------------------------------------------------------------------------\n');
fprintf('Importing Velocity Lookup Table Reference Data... \n');

fprintf('--\n');
VelLookupTablePreference = input('Would you like to use the April 2024 Brothers Test Launch Data set as reference?\n-> (y/n): ',"s");
switch VelLookupTablePreference
    case {"y","Y",""}
        % Import data from April Brothers 2024 Launch
        AprilBrothers1 = readtable("AprilBrothers_TeleMetrum_11445.csv",'VariableNamingRule','preserve'); % TeleMetrum 1
        AprilBrothers2 = readtable("AprilBrothers_TeleMetrum_11439.csv",'VariableNamingRule','preserve'); % TeleMetrum 2
        
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
    case {"n","N"}
        % prompt user for data
        LookupTableFile = input('->  Enter name of file to be used as a lookup table: ',"s");
        LookupTable1 = readtable(LookupTableFile,'VariableNamingRule','preserve');
        % !!! IMPORTANT: you will probably have to write code here to get
        % custom lookup tables!
    otherwise
        error('Yo what? I said (y/n).');
end

% Create Velocity Lookup Table
% velocity lookup table gives the PID a target velocity to achieve based
% upon a given altitude. the idea is to have the PID system control the
% velocity instead of a height. Version 4 used a target altitude of
% 10,000ft instead of a velocity lookup table

try 
    LookupRes = 0.01; % resolution of the lookup table in meters
    LookupExtend = 1.05; % extends the range of lookup table by +5%
    LookupPolyDegree = 5; % 5th degree polynomial
    
    % height
    VelLookup(:,2) = (0:LookupRes:max(LookupTable1.speed*LookupExtend))'; % create height vector with LookupRes resolution
    % velocity
    VelLookup(:,1) = polyval(...
                            polyfit(...                 % find polynomial between
                                LookupTable1.speed,...  % speed
                                LookupTable1.height,... % height
                                LookupPolyDegree...     % degree
                                   ),...
                            VelLookup(:,2)... % evaluate polynomial at desired height values
                            );
catch
    error('Error likely with importing custom data set. You likely need to write some code yourself to import a custom profile.');
end

fprintf('\nImporting Velocity Lookup Table Reference Data Complete!\n');

%% Additional Script Inputs -----------------------------------------------

fprintf('\n---------------------------------------------------------------------------\n');
fprintf('Getting Additional Script Inputs...\n');

% Determine coast phase begin time and index number
BurnoutTime = RocketEvent.Time(4); % define time of burnout
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time

% BEAVS Characteristic Parameters

fprintf('For the following parameters, press "ENTER" to select default values.\n');

fprintf('--\n');
Blade.Wdth = input('->  Blade width (meters): ');
if isempty(Blade.Wdth)
    Blade.Wdth = 7*10^-2; % convert to meters
    fprintf('Default blade width selected: %0.3f meters\n',Blade.Wdth);
end

fprintf('--\n');
Blade.Cnt = input('->  Blade count: ');
if isempty(Blade.Cnt)
    Blade.Cnt = 2; % 2 blades
    fprintf('Default blade count selected: %0.0f blades\n',Blade.Cnt);
end

fprintf('--\n');
Blade.ExtnRate = input('->  Blade extension rate (m/s): ');
if isempty(Blade.ExtnRate)
    Blade.ExtnRate = 110e-3; % meters per second
    fprintf('Default extension rate selected: %0.3f m/s\n',Blade.ExtnRate);
end

% possible pysical hardstop options are at
%   20, 25, 40, 50, 75, 60, and 100 percent
% of 60mm
fprintf('--\n');
HardstopLoc = input('Do you have a hard stop?\n If so, what % of max extension is it located at?\n If not, press "ENTER" to skip.\n0-100 corrolating to a percent\n-> ');
if isempty(HardstopLoc)
    HardstopLoc = 100; % 100%
    fprintf('Default hard stop location selected: %0.0f%%\n',HardstopLoc);
end
while (HardstopLoc <= 0) || (HardstopLoc > 100)
    HardstopLoc = input('Hardstop must be between 0 and 100\n->  ');
end

fprintf('--\n');
Blade.ExtnMAX = input('->  Blade max extension (meters): ');
if isempty(Blade.ExtnMAX)
    Blade.ExtnMAX = 0.06; % max extension of BEAVS in m
    fprintf('Default max extension selected: %0.3f meters\n',Blade.ExtnMAX);
end

% account for the hardstop
Blade.ExtnMAX = Blade.ExtnMAX*(HardstopLoc/100);

% Create Blade.Extn Vector
% used for maximum braking scenario
Blade.Extn = zeros(numel(RocketData.Time),1); % extension of BEAVS over TIME (convert to meters)
Blade.Extn(1:iterStart) = 0; % no extension in boost phase
% ramp up BEAVS extension until fully extended at extension rate
for i = iterStart:numel(Blade.Extn)-1
    % if max extension is hit set to max
    if Blade.Extn(i) > Blade.ExtnMAX
        Blade.Extn(i) = Blade.ExtnMAX;
    end
    % Extend blades
    Blade.Extn(i+1) = Blade.Extn(i) + Blade.ExtnRate*( RocketData.Time(i+1) - RocketData.Time(i) );
end

% NOISE
% see Noise.m
fprintf('--\n');
NoiseDeg = input('Degree of noise to use in simulations\n0-100 corrolating to a percent\n-> ');
if isempty(NoiseDeg)
    NoiseDeg = 0;
    fprintf('Default noise selected: %0.0f\n',NoiseDeg);
end
while (NoiseDeg < 0) || (NoiseDeg > 100)
    NoiseDeg = input('Noise must be between 0 and 100\n-> ');
end

% target altitude
fprintf('--\n');
hTarg = input('->  Target altitude (ft):');
if isempty(hTarg)
    hTarg = 10000*0.3048; % default 10,000 ft
    fprintf('Default target altitude selected: %0.0f ft\n',hTarg/0.3048);
end

fprintf('\nGetting Additional Script Inputs Complete!\n');

%% SIMULATIONS ------------------------------------------------------------
% Use Forward Euler to calculate velocity and altitude

tic
fprintf('\n---------------------------------------------------------------------------\n');
fprintf('Simulating Scenarios...\n');

% various constants
Cd.openRocket = RocketData.DragCoefficient; % OpenRocket total Cd
Cd.sim = zeros(size(Cd.openRocket)); % used as output destination from sims
A.ref = RocketData.ReferenceArea(1).*(0.01).^2; % convert to m2
rho = (RocketData.AirPressure.*100)./287.058./(RocketData.AirTemperature + 273.15); % convert to kg/m3

% Pure OpenRocket Simulation, No BEAVS
SimNum = 1;
SimName(SimNum) = {'OpenRocket'};
Time = RocketData.Time;
h    = RocketData.Altitude;
V    = RocketData.VerticalVelocity;

% Const. Braking
SimNum = 2;
SimName(SimNum) = {'Const. Braking'};
[Time(:,SimNum), h(:,SimNum), V(:,SimNum), Cd.sim(:,SimNum), Fb(:,SimNum), Fn(:,SimNum)] = ...
    FEuler(RocketData,RocketEvent,Cd.openRocket,Blade.Wdth,Blade.Cnt,Blade.ExtnRate,Blade.ExtnMAX,NoiseDeg,"No Feedback",Blade.Extn);
Fb(700:end,:) = 0;
Fn(700:end,:) = 0;
fprintf('%s Complete!\n',SimName{SimNum});

% Prediction
% !!! NOT A GREAT METHOD, UNSUPPORTED !!!
SimNum = 3;
SimName(SimNum) = {'Prediction'};
% [Time(:,SimNum), h(:,SimNum), V(:,SimNum), Cd.sim(:,SimNum), PredBladeExtn, PredBladeExtnDesire] = ...
%     FEuler(RocketData,RocketEvent,Cd.openRocket,Blade.Wdth,Blade.Cnt,Blade.ExtnRate,Blade.ExtnMAX,NoiseDeg,"Pred");
% fprintf('%s Complete!\n',SimName{SimNum});

% PID
SimNum = 4;
SimName(SimNum) = {'PID'};
Kp = 0.200e-04; % Kp = 1.200e-04; % Kp = 1.320e-04;
Ki = 1.250e-09; % Ki = 1.000e-09; % Ki = 1.000e-07;
Kd = 7.500e-05; % Kd = 6.500e-05; % Kd = 1.004e-06;
% print K values
fprintf('PID Constants:\n    Kp = %0.3e\n    Ki = %0.3e\n    Kd = %0.3e\n',Kp,Ki,Kd);
[Time(:,SimNum), h(:,SimNum), V(:,SimNum), Cd.sim(:,SimNum), PidBladeExtn, PIDu, err, Vtarg] = ...
    FEuler(RocketData,RocketEvent,Cd.openRocket,Blade.Wdth,Blade.Cnt,Blade.ExtnRate,Blade.ExtnMAX,NoiseDeg,"PID",VelLookup,Kp,Ki,Kd);
fprintf('%s Complete!\n',SimName{SimNum});

fprintf('\nSimulating Scenarios Complete!\n');
toc

%% FIGURES ----------------------------------------------------------------

fprintf('\n---------------------------------------------------------------------------\n');
fprintf('Creating Figures...\n');
tic

%% Figure 1
% Control Range Graphic

figure(fig1);
hold on
title('Control Range Graphic');
xlabel('Time (s)');

% LEFT AXIS ---------------------------------------------------------------
yyaxis left % Altitude
plot(Time(:,1:SimNum),h*3.28084);
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
plot(Time(:,1:SimNum),Cd.sim(:,1:SimNum));
ylabel('Cd total');
lims = [0 50 0 2]; axis(lims);

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','southeast');
grid on

fprintf('Fig 1 Complete\n');

%% Figure 2
% Velocity vs. Time

figure(fig2);
hold on
title('Velocity vs. Time');
xlabel('Time (s)');

% Velocity
plot(Time(:,1:SimNum),V*3.28084);
ylabel('Velocity (ft/s)');
lims = [0 50 0 1125.33]; axis(lims);

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend(SimName,'Location','SouthEast');
grid on

fprintf('Fig 2 Complete\n');

%% Figure 3
% Extension & Desired Extension

figure(fig3);
hold on
title('Extension & Desired Extension');
xlabel('Time (s)');

% prediction method is disabled
% plot(Time(:,1),[PredBladeExtn PredBladeExtnDesire]*1000);
plot(Time(:,1),PidBladeExtn*1000);
lims = [0 25 0 Blade.ExtnMAX*1000]*1.1; axis(lims);
ylabel('Extension (mm)');

yline(Blade.ExtnMAX*1000,'--','Max Extension');
legend('Actual Extn','Desired Extn','Location','SouthEast');
grid on

fprintf('Fig 3 Complete\n');

%% Figure 4
% PID Graphic

figure(fig4)
hold on
title('PID Graphic');
xlabel('Time (s)');

% plot target velocity, actualy velocity from PID data set, and error
% between the two
plot(Time(:,1),[Vtarg V(:,4) err]*3.28084);
lims = [0 40 -200 1000]; axis(lims); % cut off extraneous data
ylabel('Velocity (ft/s)');

% find time at which the target altitude was reached
TimeTrgtRched = Time(h(:,4)>=hTarg,4);
if isempty(TimeTrgtRched) % didn't get to target apogee
    TimeTrgtRched = Time(h(:,4)==max(h(:,4)),4);
    hTrgtRched = max(h(:,4));
else % passed target apogee
    TimeTrgtRched = TimeTrgtRched(1);
    hTrgtRched = hTarg;
end

% find time at which apogee was reached
TimeApogeeRched = Time(h(:,4)==max(h(:,4)),4);
hApogeeRched = max(h(:,4));

% target apogee line (black)
xline(TimeTrgtRched,'--',"Target",...
      'LabelHorizontalAlignment','left');
text(TimeTrgtRched,lims(4)-((lims(4)-lims(3))*0.1),sprintf(' t=%0.1fs\n h=%0.0fft',TimeTrgtRched,hTrgtRched*3.28084));

% actual apogee line (green)
xline(TimeApogeeRched,'--',"Apogee",...
      'LabelHorizontalAlignment','left',...
      'Color',"#77AC30",...
      "LabelVerticalAlignment","bottom");
text(TimeApogeeRched,lims(3)+((lims(4)-lims(3))*0.1),sprintf(' t=%0.1fs\n h=%0.0fft',TimeApogeeRched,hApogeeRched*3.28084),...
     'Color',"#77AC30");

legend('Target','Actual','Error','Location','NorthEast');
grid on

fprintf('Fig 4 Complete\n');

%% Figure 5
% Drag Seperation Forces

figure(fig5);
hold on;

title('Drag Seperation Forces');
xlabel('Time (s)');

plot(Time(:,1),abs(Fb(:,2))); % Force on singular BEAVS blade
plot(Time(:,1),abs(Fn(:,2))); % Force on coupler
lims = [0 30 0 2000]; axis(lims);
ylabel('Force (N)');

yyaxis right
plot(Time(:,1),Blade.Extn*1000);
ylabel('Blade Extension (mm)');

% Motor Cutoff and Legend
xline(BurnoutTime,'-',{'Motor Cutoff'},'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
legend({'Blade Force', 'Coupler Force','Blade Extension'},'Location','NorthEast');
grid on

Blade.ExtnALLOW = InterpCdInverse(2.*abs(Fb(:,2))./rho./A.ref./(V(:,2).^2),A.ref)./Blade.Wdth.*1000; % convert to mm

fprintf('Fig 5 Complete\n');

%% Figure 6
% Compare to Real
% Similar to Control Range Graphic however, with data from April Brothers
% test launch plotted as well

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

grid on;

% LEFT AXIS -----------------------------------------------------------
yyaxis right % Velocity
plot(Time, V*3.281);
plot(AprilBrothers1.time, AprilBrothers1.speed*3.281,'-','LineWidth',2);
plot(AprilBrothers2.time, AprilBrothers2.speed*3.281,'--','LineWidth',2);

ylabel('Velocity (ft/s)');
lims = [0 50 -1000 1500]; axis(lims);

legend([SimName, 'April Brothers 1', 'April Brothers 2'],'Location','SouthEast');

fprintf('Fig 6 Complete\n');

%% Figure 7

figure(fig7); hold on; grid on;

plot(VelLookup(:,2)*3.28084,    VelLookup(:,1)*3.28084         ); % plot polynomial
plot(LookupTable1.speed*3.28084,LookupTable1.height*3.28084    ); % plot data set 1
plot(LookupTable2.speed*3.28084,LookupTable2.height*3.28084,':'); % plot data set 2
title('Altitude vs. Velocity Lookup Table');
xlabel('Velocity (ft/s)'); ylabel('Altitude (ft)');
legend('Polynomial Lookup','April Brother Data Set 1','April Brother Data Set 2');
lims = axis; lims(1) = 0; axis(lims);

fprintf('Fig 7 Complete\n');

%% Conclude Figures

figure(fig1);
fprintf('\nCreating Figures Complete!\n');
toc

%% Summary of Simulation

% output summary to cmd window for each simulation
for i = 1:SimNum
    % comment out to reduce cmd window spam
    %Summary(Time(:,i), h(:,i), V(:,i), Cd.sim(:,i), RocketData, RocketEvent, SimName{:,i})
end

%% -------------------- !!! TERMINATE SCRIPT !!! --------------------
% comment out or delete to run the rest
return

%% Drag Seperation

% OPENROCKET CD's
Cd.nose   = 0.554;
Cd.nose_s = 0.034;
Cd.fore   = 0.132;
Cd.nac    = 0.002;
Cd.aft    = 0.273;
Cd.strake = 0.011;
Cd.fin    = 0.012;

% Assemble Cd's
Cd.fore_tot = Cd.nose + Cd.nose_s + Cd.fore + Cd.nac;
Cd.aft_tot = Cd.aft + Cd.strake + Cd.fin;

A.BEAVS = Blade.ExtnMAX*Blade.Wdth*Blade.Cnt;
Cd.BEAVS = InterpCd(A.BEAVS,A.ref);

% Forces
F_fore = 0.5*rho(iterStart)*V(iterStart,1).^2*Cd.fore_tot*A.ref;
F_aft  = 0.5*rho(iterStart)*V(iterStart,1).^2*(Cd.aft_tot + 2*Cd.BEAVS*(A.BEAVS/A.ref))*A.ref;

% If positive, AFT drag is LARGER than FORE (BAD!!!)
% If negative, AFT drag is LARGER than FORE (GOOD!!!)
fprintf('\n---------- Drag Seperation Forces ----------\n');
fprintf('Max BEAVS Extension: %0.0f%%\n',HardstopLoc);
fprintf('Fore Drag: %0.2f N\n',F_fore);
fprintf('Aft Drag: %0.2f N\n',F_aft);
fprintf('If positive, AFT drag is LARGER than FORE (BAD!!!)\n');
fprintf('If negative, AFT drag is LARGER than FORE (GOOD!!!)\n');
fprintf('Difference: %0.2f N\n',F_aft-F_fore);

