%% RKutta()
function [Time, h, V, Cd] = RKutta(RocketData, RocketEvent, Extension, Cd_rocket, Cd_BEAVS, BladeWdth, BladeCnt)
%% !!! WARNING THIS FUNCTION IS NO LONGER USED AN DEPRECIATED !!!

error('Function RKutta() is depreciatred and should not be used.');

% This function takes OpenRocket Data, BEAVS blade extension, and Cd's
  % and applies ode45() aka Runge-Kutta (4,5) to calculate altitude and
  % velocity over time
% INPUTS ------------------------------------------------------------------
    % RocketData   | Table from OpenRocket
    % RocketEvent  | Table of events from OpenRocket
    % Extension    | vector containing BEAVS blade extension over the
    %                OpenRocket Time vector
    % Cd_rocket    | desired coefficient of drag for the rocket w/o BEAVS
    %                can be vector w/ respect to OpenRocket Time
    % Cd_BEAVS     | desired coefficient of drag for BEAVS baldes ALONE
    %                can be vector w/ respect to OpenRocket Time
% OUTPUTS -----------------------------------------------------------------
    % Time | time vector, the same as OpenRocket time vector, just helpful
    % h    | altitude vector w/ respect to Time
    % V    | velocity vector w/ respect to Time
    % Cd   | calculated Cd vector for reference
    
    % Reference Area
    A.RefOpen = RocketData.ReferenceArea.*(0.01).^2; % convert to m

    % Calculate BEAVS Area
    A.BEAVS = Extension.*BladeWdth.*BladeCnt; % calculate BEAVS area

    % Unused
%   A.Fins = 0.00508*0.10795;
%   A.Tube = 0.02043;
%   A.Nacelle = (1.25*2.54)*(1/100)^2;
%   A.RefBEAVS = A.Tube + A.BEAVS;

    % Calculate Cd TOTAL
%   Cd = RocketData.DragCoefficient; Cd(isnan(Cd)) = 0;
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
    
    % Constants
    g = RocketData.GravitationalAcceleration;
    Rd = 287.058;
    rho = (RocketData.AirPressure.*100)./Rd./(RocketData.AirTemperature + 273.15);
    m = RocketData.Mass/1000;
%   T = RocketData.Thrust;
    
    tspan = Time(iterStart:end); % Solution output at TIME
    h0 = [h(iterStart) V(iterStart)]; % initial conditions
    [~, H] = ode45(@(t,h) odesys(t,h,A.RefOpen,Cd,g,rho,m,Time,FeedbackBool), tspan, h0);
    
    h(iterStart:end) = H(:,1);
    V(iterStart:end) = H(:,2);
end

% ODE System for ODE45
function dhdt = odesys(t,h,A,Cd,g,rho,m,Time,FeedbackBool)
    if FeedbackBool == true
%        [A, Cd] = GetFeedback([Time(i) Time(i+1)], h(i), V(i));
    else % interpolate for parameters at time 't'
        A   = interp1(Time,A  ,t);
        Cd  = interp1(Time,Cd ,t);
    end
    
    % interpolate parameters at time 't'
    g   = interp1(Time,g  ,t);
    rho = interp1(Time,rho,t);
    m   = interp1(Time,m  ,t);
    
    dhdt = [h(2);...
            -g - 1/2*rho*h(2)*abs(h(2))*Cd*A/m];
end