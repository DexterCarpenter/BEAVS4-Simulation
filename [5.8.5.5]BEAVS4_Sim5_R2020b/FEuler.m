%% FEuler()
function [Time, h, V, Cd, varargout] = FEuler(RocketData,RocketEvent,Cd_rocket,BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX,NoiseDeg,Method,varargin)
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

% For Solving SolveFdrag()
ma = 12.711; % kg
mb =  7.517; % kg
% allocate
Fb = zeros(numel(Time),1);
Fn = zeros(numel(Time),1);

%% NO FEEDBACK ------------------------------------------------------------
% Use given blade extension
% FEuler(...,"No Feedback",BladeExtn)

if Method == "No Feedback"
    % Resolve varargin
    % {1} is Blade Extension
    BladeExtn = varargin{1};
    
    % Define Variables
    % Calculate BEAVS Blade Area
    A_BEAVS = BladeExtn.*BladeWdth.*BladeCnt; % calculate BEAVS area
    
    % get BEAVS Cd
    Cd_BEAVS = InterpCd(A_BEAVS,A_ref);
    
    % get Total Cd
    Cd = Cd_rocket + Cd_BEAVS.*(A_BEAVS./A_ref);
    
    % Iterate Forward Euler
    for i = iterStart:numel(Time)-1   
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i), NoiseDeg);
        [Fb(i),Fn(i)] = SolveFdrag(ma,mb,RocketData.TotalAcceleration(i),g(i),RocketData.DragForce(i));
    end
    
    varargout{1} = Fb;
    varargout{2} = Fn;
    
    return
end


%% PREDICTION -------------------------------------------------------------
% predict apogee and get a desired blade extension
% [...,BladeExtn,BladeExtnDesire] = FEuler(...,"Pred", N/A )

if Method == "Pred" || Method == "Predict"
    % Resolve varargin
    
    % no extra args needed
    
    % Define Variables
    % Allocate BladeExtn and BladeExtnDesire vector
    BladeExtn       = zeros(numel(Time),1);
    BladeExtnDesire = zeros(numel(Time),1);
    Cd = Cd_rocket;
    
    % Iterate Forward Euler
    for i = iterStart:numel(Time)-1
        % Get Cd and Blade Extension
        [Cd(i+1),BladeExtn(i+1),BladeExtnDesire(i+1)] = ...
            GetFeedbackPred([Time(i) Time(i+1)],h(i),V(i),g(i),rho(i),A_ref(i),m(i),a(i),BladeExtn(i),Cd_rocket(i),BladeWdth,BladeCnt,BladeExtnRate,BladeExtnMAX);
        % Step FEuler
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i), NoiseDeg);
    end
    
    % Resolve varargout
    varargout{1} = BladeExtn;
    varargout{2} = BladeExtnDesire;
    
    return
end


%% PID Feedback -----------------------------------------------------------
% use PID loop to control blade extension
% FEuler(...,"PID",VelLookup,Kp,Ki,Kd)

if Method == "PID"
    % Resolve varargin
    VelLookup = varargin{1};
    Kp = varargin{2};
    Ki = varargin{3};
    Kd = varargin{4};
    
    % Define Variables
    BladeExtn = zeros(numel(Time),1);
    Cd        = Cd_rocket;
    u         = zeros(numel(Time),1);
    err       = zeros(numel(Time),1);
    Vtarg     = zeros(numel(Time),1);
    
    % Iterate Forward Euler
    for i = iterStart:numel(Time)-1
        
        % Get new Target
        Vtarg(i) = interp1(VelLookup(:,1), VelLookup(:,2), h(i));
        %Vtarg(i) = 1; % overwrite to step response
        if isnan(Vtarg(i)) == false
            % Update PID Control Function
            [u(i), err(i)] = GetFeedbackPID([Time(i) Time(i+1)], Kp, Ki, Kd, u(i-1), [V(i-2) V(i-1) V(i)], Vtarg(i));
        else
            % Else set control function to zero (no extension)
            u(i) = 0;
        end
        
        % Adjust Blade Extension
        [Cd(i+1),BladeExtn(i+1)] = UpdateExtension([Time(i) Time(i+1)], A_ref(i), Cd_rocket(i), BladeExtn(i), BladeWdth, BladeCnt, BladeExtnRate, BladeExtnMAX, u(i));
        
        % Step FEuler
        [h(i+1), V(i+1)] = StepFEuler([Time(i) Time(i+1)], h(i), V(i), g(i), rho(i), Cd(i), A_ref(i), m(i), NoiseDeg);
        
        % Calculate forces
%         [Fb(i),Fn(i)] = SolveFdrag(ma,mb,RocketData.TotalAcceleration(i),g(i),RocketData.DragForce(i));
    end
    
    % Resolve varargout
    varargout{1} = BladeExtn;
    varargout{2} = u;
    varargout{3} = err;
    varargout{4} = Vtarg;
%     varargout{3} = Fb;
%     varargout{4} = Fn;
    
    return
end

%% Method not Recognized --------------------------------------------------

fprintf('\n ERROR in FEuler\n\nMethod "%s" is not recognized by FEuler\n',Method);

end

%% StepEuler()
function [hout, Vout] = StepFEuler(Time, h, V, g, rho, Cd, A_ref, m, NoiseDeg)
% Calculate V
Vout = V + (-g - 1/2*rho*V*abs(V)*Cd*A_ref/m)*(Time(2) - Time(1));
% Calculate h
hout = h + V*(Time(2) - Time(1));

Vout = Noise([V Vout],NoiseDeg);
hout = Noise([h hout],NoiseDeg);
end



