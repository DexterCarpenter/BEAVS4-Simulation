%% FEuler()
function [Time, h, V, Cd, varargout] = FEuler(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX,Method,varargin)
%% Universal
% the following variables are used by each Method

% Determine where to begin iteration (start at burnout)
BurnoutTime = RocketEvent.Time(4);
iterStart = find(RocketData.Time==BurnoutTime); % iteration start time

% Allocate Vectors
Time           = RocketData.Time;                          % mirror
h              = zeros(numel(Time),1);                     % total vector
h(1:iterStart) = RocketData.Altitude(1:iterStart);         % boost phase
V              = zeros(numel(Time),1);                     % total vector
V(1:iterStart) = RocketData.VerticalVelocity(1:iterStart); % boost phase

% Constants
g = RocketData.GravitationalAcceleration; % m/s2
rho = (RocketData.AirPressure.*100)./287.058./(RocketData.AirTemperature + 273.15); % convert to kg/m3
m = RocketData.Mass/1000; % convert to kg
A_ref = RocketData.ReferenceArea.*(0.01).^2; % convert to m2
a = deg2rad(RocketData.VerticalOrientation); % NOT ANGLE OF ATTACK!!! (rad)

%% No Feedback ------------------------------------------------------------
% Use given blade extension
% FEuler(...,"No Feedback",Extn)

if Method == "No Feedback"    
    %% Resolve varargin
    % {1} is Blade Extension
    BladeExtn = varargin{1};
    
    %% Define Variables
    % Calculate BEAVS Blade Area
    A_BEAVS = BladeExtn.*BladeWdth.*BladeCnt; % calculate BEAVS area
    
    % get BEAVS Cd
    Cd_BEAVS = InterpCd(A_BEAVS,A_ref);
    
    % get Total Cd
    Cd = Cd_rocket + Cd_BEAVS.*(A_BEAVS./A_ref);
    
    %% Iterate Forward Euler
    for i = iterStart:numel(Time)-1   
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i));
    end
    return
end


%% Prediction -------------------------------------------------------------
if Method == "Pred"
    %% Resolve varargin
    
    
    %% Define Variables
    % Allocate Extn and ExtnDesire vector
    BladeExtn  = zeros(numel(Time),1);
    ExtnDesire = zeros(numel(Time),1);
    Cd = Cd_rocket;
    
    %% Iterate Forward Euler
    for i = iterStart:numel(Time)-1   
        % Get Cd and Blade Extension
        [Cd(i+1),BladeExtn(i+1),ExtnDesire(i+1)] = GetFeedbackPred([Time(i) Time(i+1)],h(i),V(i),g(i),rho(i),A.RefOpen(i),m(i),a(i),Extn(i),Cd_rocket(i),BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX);
        % Step FEuler
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i));
    end
    return
end



%% PID Feedback -----------------------------------------------------------
if Method == "Pred"
    %% Resolve varargin
    
    
    %% Define Variables
    
    
    %% Iterate Forward Euler
    for i = iterStart:numel(Time)-1   
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i));
    end
    return
end




end

%% StepEuler()
function [h, V] = StepFEuler(Time, h, V, g, rho, Cd, A_ref, m)
% Calculate V
V = V + (-g - 1/2*rho*V*abs(V)*Cd*A_ref/m)*(Time(2) - Time(1));
% Calculate h
h = h + V*(Time(2) - Time(1));
end



