function [Time, h, V, Cd, Extn, ExtnDesire] = FEulerFeedback(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BEAVSExtnMAX,varargin)
%% !!! WARNING THIS FUNCTION IS NO LONGER USED AN DEPRECIATED !!!

error('Function FEulerFeedback() is depreciatred and should not be used.');

%% FEulerFeedback()
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
    % FeedbackBool | true - use GetFeedback.m
    %              | false - do not use GetFeedback.m
% OUTPUTS -----------------------------------------------------------------
    % Time | time vector, the same as OpenRocket time vector, just helpful
    % h    | altitude vector w/ respect to Time
    % V    | velocity vector w/ respect to Time
    % Cd   | calculated Cd vector for reference
    
% Determine whether to use PID for feedback
if isempty(varargin) == false
    UsePID = true;
    K = varargin{1};
else
    UsePID = false;
end

% Rocket Reference Area
A.RefOpen = RocketData.ReferenceArea.*(0.01).^2; % convert to m

% Determine where to begin iteration (start at burnout)
BurnoutTime = RocketEvent.Time(4);
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time

% Allocate vectors for Beavs Data
Time           = RocketData.Time;
h              = zeros(numel(Time),1);             % total vector
h(1:iterStart) = RocketData.Altitude(1:iterStart); % boost phase
V              = zeros(numel(Time),1);                     % total vector
V(1:iterStart) = RocketData.VerticalVelocity(1:iterStart); % boost phase
a              = deg2rad(RocketData.VerticalOrientation); % NOT ANGLE OF ATTACK!!!

% Constants
g = RocketData.GravitationalAcceleration;
Rd = 287.058;
rho = (RocketData.AirPressure.*100)./Rd./(RocketData.AirTemperature + 273.15);
m = RocketData.Mass/1000;
T = RocketData.Thrust;

% Allocate Extn and ExtnDesire vector
Extn       = zeros(numel(Time),1);
ExtnDesire = zeros(numel(Time),1);
Cd = Cd_rocket;

for i = iterStart:numel(Time)-1
    if UsePID == false
        [Cd(i+1),Extn(i+1),ExtnDesire(i+1)] = GetFeedbackPred([Time(i) Time(i+1)],h(i),V(i),g(i),rho(i),A.RefOpen(i),m(i),a(i),Extn(i),Cd_rocket(i),BladeWdth,BladeCnt,BladeExtnRate,BEAVSExtnMAX);
    else
        [Cd(i+1),Extn(i+1),ExtnDesire(i+1)] = GetFeedbackPID([Time(i) Time(i+1)],h(i),V(i),g(i),rho(i),A.RefOpen(i),m(i),a(i),Extn(i),Cd_rocket(i),BladeWdth,BladeCnt,BladeExtnRate,BEAVSExtnMAX,K);
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