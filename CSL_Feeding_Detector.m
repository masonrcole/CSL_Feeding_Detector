%CSL_Feeding_Detector
%
%   Detects prey capture by California sea lions using the surge- and
%   heave-axis acceleration (from 3-axis acceleration data), as validated
%   by Cole et al. (in Review)
%
%   INPUT PARAMETERS:
%       A: (required) 3-axis acceleration data from a head-mounted
%           accelerometer, organized in 3 columns as [sway surge heave].
%           Be sure acceleration data are in units 'g' (1g = 9.81 m
%           s^-2). Surge-axis acceleration data must be sampled (or
%           arranged) such that positive values indicate forward (direction
%           of swimming) acceleration, and negative values indicate
%           backward acceleration. 
%       rateS: (required) tag sampling rate, in Hz
%       rateA: (optional) desired analysis rate (Hz), if different than
%           rateS. 
%       time: (optional) timestamps corresponding to A. Required IF
%           plotting prey capture against depth or if using a depth threshold
%       depth: (OPTIONAL) depth vector; add if a plot of feeding events
%           overlain on a depth profile is desired
%       rateD: (OPTIONAL; Mandatory if plotting depth) depth sampling rate.
%       depth_time: (OPTIONAL) timestamps corresponding to optional depth
%           vector; essential if depth sampling rate differs from that of
%           acceleration data
%       depth_threshold: (OPTIONAL) depth in meters (positive values); if
%           depth_threshold is used, only detected feeding events
%           deeper than the threshold will appear in plots and outputs. 
%
%   OUTPUTS:
%       Summary: summary statistics describing each detected feeding
%           event. Arranged as:
%               Column 1: Event #
%               Column 2: Detected event Duration
%               Column 3: Maximum value of smoothed heave-axis Jerk
%               Column 4: Maximum (negative) value of dynamic surge-axis deceleration
%               Column 5: Integral of smoothed heave-axis Jerk (area under
%                        heave-Jerk curve in detected zone)
%               Column 6: Integral of |dynamic surge-axis values| in
%                         detected zone
%               Column 7: timestamp of onset of event
%
%       Event_Data: relevent data for each event, all in one big matrix
%               Column 1: Feeding Event #
%               Column 2: Feeding Event Data timestamps
%               Column 3: smoothed heave-axis Jerk data
%               Column 4: filtered (~dynamic) surge-axis acceleration data
%
%       Plots: Fig. 1: Surge-axis filtered acceleration (top) and
%                       Heave-axis smoothed Jerk (bottom), for visual
%                       inspection of acceleration signals
%              Fig. 2: (if depth is added as an input): detected prey capture events plotted on depth profile 
%

%------------------------------------------------------------------------------------------------------------
% DEFINE PARAMETERS HERE:
%    It's important that the variable names do not change, as they are
%    inputed into a function

A = [acc(:,1) acc(:,2) acc(:,3)];  % [sway surge heave]
rateS = 50;
rateA = 50;                             % if rateA not specified, set rateA = [];
time = t;         % if time not specified, set time = [];
depth = depth;               % if depth not added, set depth = [];
rateD = 1;                             % depth sampling rate in Hz. if rateD not specified, set rateD = [];
depth_time = ptime;    % if depth_time not specified, set depth_time = [];
depth_threshold = 5;                    % if depth_threshold not specified, set = [];

% Now Run Script!
%-------------------------------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------------------
%internal function: CSL_Feeding
[Summary,Event_Data] = CSL_Feeding(A,rateS,rateA,time,depth,rateD,depth_time,depth_threshold);
