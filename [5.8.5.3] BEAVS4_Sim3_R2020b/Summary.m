%% Summary()
function Summary(Time, h, V, Cd, RocketData, RocketEvent, Name)
% Format and Summarize Data to Command Window
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