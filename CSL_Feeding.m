function [Summary,Event_Data] = CSL_Feeding(A,rateS,rateA,time,depth,rateD,depth_time,depth_threshold) 
% NEW CSL_Feeding_ID
%
%   Identifies prey capture attempts from head-mounted 2- or 3-axis
%   accelerometers (heave and surge) in California sea lions. 
%
%
%   [note from Mason to Mason: MUST INCLUDE movingmean.m IN PACKAGE]
%
%   INPUTS: 
%       A: 3-axis acceleration; axes must be vectors of equal length,
%           arranged as [sway surge heave]. Sway is not used, so a column
%           of place-holder values may be used in its place if desired
%       rateS: sampling rate of the tag (Hz)
%       rateA: (OPTIONAL; Mandatory if plotting against depth) analysis rate for prey capture analyses (Hz), if less than rateS. If rateA
%           is not specified, default rateA = rateS
%       time: (OPTIONAL; Mandatory if plotting against depth) timestamps corresponding to acceleration data;
%           vector must be same length as acceleration vector. Useful in
%           matching feeding events to depth (if depth sampled slower,
%           etc.)
%       depth: (OPTIONAL) depth vector; add if a plot of feeding events
%           overlain on a depth profile is desired
%       rateD: (OPTIONAL; Mandatory if plotting depth) depth sampling rate.
%       depth_time: (OPTIONAL) timestamps corresponding to optional depth
%           vector; essential if depth sampling rate differs from
%           acceleration data
%       depth_threshold: (OPTIONAL) depth in meters (positive values); if
%           depth_threshold is used, only detected feeding events
%           deeper than the threshold will appear in plots and outputs.
%           
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
%               Column 3: smoothed heave-Jerk data
%               Column 4: dynamic surge acceleration data
%
%       Plots: Fig. 1: Surge-axis filtered acceleration (top) and
%                       Heave-axis smoothed Jerk (bottom), for visual
%                       inspection of acceleration signals
%              Fig. 2: (if depth is added as an input): detected prey capture events plotted on depth profile 
%
%----------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
% Thresholds that are independent of sampling rate:
YDYNthreshold = -0.7;                    % surge-axis deceleration threshold
surgeTH = 1.0;                         % surge-axis acceleration threshold

%-----------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------
% (Internal processes):


% Calculate acceleration matrix, jerk, and time variables for chosen Analysis Rate

ARin = isempty(rateA);
if ARin == 0                    % if rateA is not empty, rateA = fs
    fs = rateA;
else
    rateA = rateS;              % if rateA is empty, rateA = rateS = fs
    fs = rateA;
end
fac = rateS/rateA;
A = A(1:fac:size(A,1),:);       % subsamples raw acceleration data (sampling rate) into specified Analysis rate
j = (9.81*fs)*sqrt(diff(A).^2*ones(3,1));       % j = calculated jerk from A with sampling speed fs
t_jerk = [0:1/(fs):length(j)]';                 %t_jerk calculates the (time) x-axis for the jerk figures
t_jerk = t_jerk(1:length(j),1);

%------------------------------------------------------------------------------------------------------------------------
% create Accel timestamps as 4th column of A (so A(:,4) is based of SR and AR)
ratestep = 1/fs;
for loop = 1:size(A,1)
    if loop == 1
        A(loop,4) = 0;
    else
        A(loop,4) = (loop.*ratestep)-ratestep;
    end
end

timecheck = isempty(time);
% if timecheck == 1
% %     time = time(1:fac:size(time,1),:);              %Doesn't make sense? IDEA: make A(:,4) same whether 'time' is included or not. BUT use 'time' for syncing with depth
%     A(:,4) = time;
% end

%------------------------------------------------------------------------------------------------------------------------------
% if depth is an input & depth_time is not, create a depth_time vector 'D_T'.
%   Needed if a depth_threshold is imposed
Dchk = isempty(depth);
DTchk = isempty(depth_time);
DThresh_check = isempty(depth_threshold);
if Dchk == 0 & DTchk == 1 & DThresh_check == 0
    rateD_chk = isempty(rateD);
    if rateD_chk == 0
        %if rateD given, use this (compared with rateA) to create the vector
        if timecheck == 0
            D_T = (time(1):(1/rateD):time(end));
            D_T = D_T(1:length(depth));
        else
            for lp = 1:size(depth,1)
                if lp == 1
                    D_T(lp) = 0;
                else
                    D_T(lp) = lp.*(1/rateD) - rateD;
                end
            end
        end
    else
        Lfac = size(A,1)/size(depth,1);         %should be the difference between depth and accel sampling rates
        rD = rateA/Lfac;
        if timecheck == 0
            D_T = (time(1):(1/rD):time(end));
            D_T = D_T(1:length(depth));        
        else
            for lp = 1:size(depth,1)
                if lp == 1
                    D_T(lp) = 0;
                else
                    D_T(lp) = lp.*(1/rD) - rD;
                end
            end
        end
    end
end
if Dchk == 0 & DTchk == 0 & DThresh_check == 0
    D_T = depth_time;
end
%--------------------------------------------------------------------------------------------------------------------------------
% Calculate heave-axis Jerk threshold ('ZJthreshold') from sampling rate using Cali's empirical data
ZJthreshold = 22.711.*((fs)^0.4327);                      % ZJthrehsold based on analysis (or sampling) rate
%------------------------------------------------------------------------------------------------------------------------
% Calculate heave jerk, smoothed heave jerk, and dynamic surge axis acceleration 

zj = (9.81*fs)*sqrt(diff(A(:,3)).^2) ;                              % zj = heave axis jerk
smZJwindow = round(fs/20);                                      % moving mean window = SampleRate/20
smoothZJ = movingmean(zj,smZJwindow);                           % smoothZJ is a moving mean (window = smZJwindow) of heave-jerk

smYDYNwindow = round(fs/2);                                     % fs/2 = moving mean window to calculate acceleration due to orientation along surge axis (to be subtracted from raw)
smoothYDYN = movingmean(A(:,2),smYDYNwindow);                   % smoothYDYN = surge acceleration due to tag's orientation
YDYN = A(:,2)-smoothYDYN;                                       % YDYN = (raw surge accel)-(orientation accel) = surge dynamic acceleration
%-----------------------------------------------------------------------------------------------------------------------
% Combine A and Jerk variables into one matrix with more columns (all the same length): named AJdata
jdata = [j(1);j];
smoothzjdata = [smoothZJ(1);smoothZJ];
t_jerkdata = [t_jerk;t_jerk(end) + ratestep];
AJdata = [A jdata smoothzjdata YDYN]; % columns: 1-3=x,y,z accel; 4=time; 5=j; 6=smoothZJ; 7=YDYN.
%------------------------------------------------------------------------------------------------------------------------
% Find the values and time stamps of smoothZJ and YDYN that surpasses calculated thresholds

ZJyesIND = find(smoothZJ > ZJthreshold);                        % index of smoothZJ values that surpass ZJ threshold
ZJyesVALUES = smoothZJ(ZJyesIND);                               % smoothZJ values that correspond to the index 
ZJyesTIME = A(ZJyesIND,4);                                      % time stamp corresponding to index

YDYNyesIND = find(YDYN < YDYNthreshold);                        % index of YDYN values that surpass YDYN threshold
YDYNyesVALUES = YDYN(YDYNyesIND);                               % corresponding values
YDYNyesTIME = A(YDYNyesIND,4);                                  % corresponding time stamp
%--------------------------------------------------------------------------------------------------------------------------
% get an index and times of all possible event times

RawTIMESind = find(smoothzjdata > ZJthreshold | YDYN < YDYNthreshold);      % index of all rows that contain EITHER smoothZJ>threshold, OR YDYN>YDYNthreshold, or both
RawTIMES = AJdata(RawTIMESind,4);                                           % corresponding timestamps
%-------------------------------------------------------------------------------------------------------------------------
% create ZJ EVENT NUMBER labels: within 1 sec of the previous point gets the same event label

ZJevents = [ZJyesTIME ZJyesTIME];                                           % create 'placeholder' matrix of correct size, 2nd column = timestamps
ZJevents(:,1)=0;                                                            % set all 1st column = 0
ZJevents(1,1)=1;                                                            % set 1st ZJyes to 1 (first event)
for timeloop = 2:size(ZJevents,1)                                           % set all events within 1 sec of the previous timepoint as being part of that event; otherwise NEXT event 
    if ZJevents(timeloop,2) < ZJevents(timeloop-1,2)+1
        ZJevents(timeloop,1) = ZJevents(timeloop-1,1);
    else ZJevents(timeloop,1) = ZJevents(timeloop-1,1)+1;
    end
end
%------------------------------------------------------------------------------------------------------------------------
% make a summary matrix for start and end times of each ZJ event

for numEvents = 1:ZJevents(end,1),                                  % make a summary matrix for start and end times of each ZJ event
    zjEventIND = find(ZJevents(:,1)==numEvents);                    % index of each seperate ZJ event
    zjStartTime = ZJevents(zjEventIND(1),2);                        %ZJ event start time
    zjEndTime = ZJevents(zjEventIND(end),2);                        %ZJ event end time
    if numEvents == 1,
        ZJeventsSummary = [numEvents zjStartTime zjEndTime];
    else
        newLine = [numEvents zjStartTime zjEndTime];
        ZJeventsSummary = [ZJeventsSummary; newLine];               % end up with matrix: c1=event #, c2=event start, c3=event end
    end
end
%-----------------------------------------------------------------------------------------------------------------------
% create YDYN EVENT NUMBER labels: within 1 sec of the previous point gets the same event label

YDYNevents = [YDYNyesTIME YDYNyesTIME];
YDYNevents(:,1)=0;
YDYNevents(1,1)=1;                                                  % set 1st ZJyes to 1 (first event)
for timeloop = 2:size(YDYNevents,1),
    if YDYNevents(timeloop,2) < YDYNevents(timeloop-1,2)+1,
        YDYNevents(timeloop,1) = YDYNevents(timeloop-1,1);
    else YDYNevents(timeloop,1) = YDYNevents(timeloop-1,1)+1;
    end
end
%-----------------------------------------------------------------------------------------------------------------------
% make a summary matrix for start and end times of each YDYN event

for numEvents = 1:YDYNevents(end,1)                                 % make a summary matrix for start and end times of each YDYN event
    ydynEventIND = find(YDYNevents(:,1)==numEvents);                % index of each seperate YDYN event
    ydynStartTime = YDYNevents(ydynEventIND(1),2);                  %YDYN event start time
    ydynEndTime = YDYNevents(ydynEventIND(end),2);                  %YDYN event end time
    if numEvents == 1
        YDYNeventsSummary = [numEvents ydynStartTime ydynEndTime];
    else
        newLine = [numEvents ydynStartTime ydynEndTime];
        YDYNeventsSummary = [YDYNeventsSummary; newLine];
    end
end
%-------------------------------------------------------------------------------------------------------------------------
% Link up events, delete YDYNeventsSummary rows (from TrueZJ) that are outside the range of ZJevents time stamps

TrueZJ = [ZJeventsSummary ones(size(ZJeventsSummary,1),1)];
INDTrueYDYN = [];
for possLOOP = 1:size(ZJeventsSummary,1)
    searchY = find(YDYNevents(:,2) > ZJeventsSummary(possLOOP,2)-0.05 & YDYNevents(:,2) < ZJeventsSummary(possLOOP,3)+0.2);  % Find time values of YDYNevents that are within ZJeventsSummary(start)-0.05sec and ZJeventsSummary(end)+0.2sec   
    if isempty(searchY)
        TrueZJ(possLOOP,4) = 0;
    end
end
removeIND = find(TrueZJ(:,4) == 0);
TrueZJ(removeIND,:) = [];
TrueZJ(:,4) = [];
%--------------------------------------------------------------------------------------------------------------------------
% Link up events, delete ZJeventsSummary rows (from TrueYDYN) that are outside the range of ZJevents time stamps

TrueYDYN = [YDYNeventsSummary ones(size(YDYNeventsSummary,1),1)];
INDTrueZJ = [];
for possLOOPy = 1:size(YDYNeventsSummary,1)
    searchZ = find(ZJevents(:,2) > YDYNeventsSummary(possLOOPy,2)-0.05 & ZJevents(:,2) < YDYNeventsSummary(possLOOPy,3)+0.2);  % Find time values of ZJevents that are within YDYNeventsSummary(start)-0.05sec and YDYNeventsSummary(end)+0.2sec   
    if isempty(searchZ)
        TrueYDYN(possLOOPy,4) = 0;
    end
end
removeINDy = find(TrueYDYN(:,4) == 0);
TrueYDYN(removeINDy,:) = [];
TrueYDYN(:,4) = [];
%--------------------------------------------------------------------------------------------------------------------------
% matching up TrueZJ and TrueYDYN

if size(TrueZJ,1) == size(TrueYDYN,1),
    TrueEvents = [TrueZJ TrueYDYN(:,2:3)];
elseif size(TrueZJ,1) > size(TrueYDYN,1),
    TrueEvents = [TrueZJ zeros(size(TrueZJ,1),2)];
    for zloop = 1:size(TrueZJ,1),
        Ymatch = find(TrueYDYN(:,2) > TrueZJ(zloop,2)-1 & TrueYDYN(:,2) < TrueZJ(zloop,2)+2); %for each ZJ start time, find the YDYN start time that is within + or - 1 seconds of that
        if ~isempty(Ymatch) & size(Ymatch,1) == 1,   % if there is 1 matching YDYN row...
            TrueEvents(zloop,4:5) = TrueYDYN(Ymatch,2:3);
        end
    end
    CUT = find(TrueEvents(:,4) == 0);
    TrueEvents(CUT,:) = [];
elseif size(TrueYDYN,1) > size(TrueZJ,1),
    TrueEvents = [TrueYDYN(:,1) zeros(size(TrueYDYN,1),2) TrueYDYN(:,2:3)];
    for yloop = 1:size(TrueYDYN,1),
        Zmatch = find(TrueZJ(:,2) > TrueYDYN(yloop,2)-1 & TrueZJ(:,2) < TrueYDYN(yloop,2)+2); %for each YDYN start time, find the ZJ start time that is within + or - 1 seconds of that
        if ~isempty(Zmatch) & size(Zmatch,1) == 1,   % if there is 1 matching ZJ row...
            TrueEvents(yloop,2:3) = TrueZJ(Zmatch,2:3);
        end
        
    end
    CUT = find(TrueEvents(:,2) == 0);
    TrueEvents(CUT,:) = [];
end
%------------------------------------------------------------------------------------------------------------------------

% Combines YDYN and ZJ event summaries into whole events (start and end times)

DurThresh = 0.05;

if ~isempty(TrueEvents),
    for i=1:size(TrueEvents,1),         % Combines YDYN and ZJ event summaries into whole events (start and end times)
        if i==1,
            FirstEvent = [TrueEvents(i,1) min(TrueEvents(i,2:5)) max(TrueEvents(i,2:5))];
            AllTrueEvents = FirstEvent;
        else
            NewEvent = [TrueEvents(i,1) min(TrueEvents(i,2:5)) max(TrueEvents(i,2:5))];
            AllTrueEvents = [AllTrueEvents; NewEvent];
        end
    end

        
%Use AllTrueEvents to: 1)build a matrix for each event, 2)identify multiple
%sub-events within an event, or remove extra data

    % columns: 1)time 2)ZJ data 3)YDYN data 

    for I = 1:size(AllTrueEvents,1),
        startIND = find(AJdata(:,4)==AllTrueEvents(I,2));
        EndIND = find(AJdata(:,4)==AllTrueEvents(I,3));
        EventData = [];
        EventData = [AJdata(startIND:EndIND,4) AJdata(startIND:EndIND,4) smoothzjdata(startIND:EndIND) YDYN(startIND:EndIND)];
        EventData(:,1) = 0;

        % from here, use indexes for each part of each event to compare values
        % for the other part 
        zjIND = find(EventData(:,3) > ZJthreshold);         % find zj data that satisfies requirements
        ydynIND = find(EventData(:,4) < YDYNthreshold);     % find ydyn data that satisfies requirements
        % all data within each 'sub-event' should be within 1)0.3sec of
        % previous qualifying data point, and 2)0.15sec of all data of
        % corresponding data type

        for D = 1:size(zjIND,1),        %zj data, proximity to previous point
            if D == 1,
                EventData(zjIND(D),1) = 1;
            elseif EventData(zjIND(D),2) < EventData(zjIND(D-1),2)+0.3, 
                EventData(zjIND(D),1) = EventData(zjIND(D-1),1);
            elseif EventData(zjIND(D),2) > EventData(zjIND(D-1),2)+0.3,
                EventData(zjIND(D),1) = EventData(zjIND(D-1),1)+1;      %next number in 1st EventData column if there is a gap
            end
        end

        EventData(:,5) = EventData(:,4);    % create a column of 1's in column 5: where sub-events in ydyn will be seperated
        EventData(:,5) = 0;
        for D = 1:size(ydynIND,1),        %ydyn data, proximity to previous point
            if D == 1,
                EventData(ydynIND(D),5) = 1;
            elseif EventData(ydynIND(D),2) < EventData(ydynIND(D-1),2)+0.3, 
                EventData(ydynIND(D),5) = EventData(ydynIND(D-1),5);
            elseif EventData(ydynIND(D),2) > EventData(ydynIND(D-1),2)+0.3,
                EventData(ydynIND(D),5) = EventData(ydynIND(D-1),5)+1;      %next number in 1st EventData column if there is a gap
            end
        end

        % we now have (in EventData): column 1 showing zj sub events, and
        % column 5 showing ydyn sub events

        % next, compare data points of each sub-event - of each type - to
        % non-zero values for the other type. For those values that don't have
        % near enough counterparts in the othe type, change sub-event value to
        % 0 (in column 1 or 5):
        for D = 1:size(zjIND,1),        %zj data, proximity to ydyn sub-event
            crossThresh = round(0.075.*fs);  %how many rows to search in both directions in ydyn subevent column
            searchmin = zjIND(D)-crossThresh;
            searchmax = zjIND(D)+crossThresh;
            if searchmin <= 0,
                searchmin = 1;
            end
            if searchmax > size(EventData,1),
                searchmax = size(EventData,1);
            end
            if 2.*crossThresh > size(EventData,1),
                searchmin = 1;
                searchmax = size(EventData,1);
            end 
            ydynsearchzone = EventData(searchmin:searchmax,5);
            zeroIND = find(ydynsearchzone ~= 0);                % find any non-zero values within window in ydyn search zone
            if isempty(zeroIND),
                EventData(zjIND(D)) = 0;
            end
        end

        %...and same for ydyn now
        for D = 1:size(ydynIND,1),        %zj data, proximity to ydyn sub-event
            crossThresh = round(0.075.*fs);  %how many rows to search in both directions in ydyn subevent column
            searchmin = ydynIND(D)-crossThresh;
            searchmax = ydynIND(D)+crossThresh;
            if searchmin <= 0,
                searchmin = 1;
            end
            if searchmax > size(EventData,1),
                searchmax = size(EventData,1);
            end 
            if 2.*crossThresh > size(EventData,1),
                searchmin = 1;
                searchmax = size(EventData,1);
            end        
            zjsearchzone = EventData(searchmin:searchmax,1);
            zeroIND = find(zjsearchzone ~= 0);                % find any non-zero values within window in ydyn search zone
            if isempty(zeroIND),
                EventData(ydynIND(D)) = 0;
            end
        end

        %NOW, we have to organize and put everything we just looped through into a variable
        %so it doesn't get lost
        EventData(:,4:6) = EventData(:,3:5);
        EventData(:,3) = EventData(:,2);
        EventData(:,2) = EventData(:,6);
        EventData(:,6) = [];
        EventData(:,2:6) = EventData(:,1:5);
        EventData(:,1) = I;
        % column 1 = Event number
        % column 2 = zj sub events
        % column 3 = ydyn sub events
        % column 4 = time
        % column 5 = zj data
        % column 6 = ydyn data
        E = EventData;

        if I == 1,
            AllData = E;
        else
            AllData = [AllData; E];
        end
    end
    
    chopIND = find(AllData(:,2) == 0 & AllData(:,3) == 0);
    AllData(chopIND,:) = [];            %remove all data that isn't included

        % define beginning of event as the beginning of Z-Jerk
    for allE = 1:size(unique(AllData(:,1)),1)
        iE = find(AllData(:,1) == allE);
        iZ = find(AllData(:,1) == allE & AllData(:,2) ~= 0);

        if ~isempty(iZ)    
            if iE(1) < iZ(1)                 %if iE starts before iZ, remove the redefine iE so that it starts when iZ starts
                iEcut = find(iE < iZ(1));
                Ecut = iE(iEcut);
                AllData(Ecut,:) = [];
            end
        else
            AllData(iE,:) = [];
        end
    end
    
    % re-number AllData(:,1) so that no event numbers are skipped (so that
    % the next part of the code works...
    list = unique(AllData(:,1));
    for L = 1:size(list,1)
        Elist = find(AllData(:,1) == list(L));
        AllData(Elist,1) = L;
    end


%-----------------------------------------------------------------------------------------------------------------------
    %Summarize events and sub-events
    m = size(unique(AllData(:,1)),1);
    AllEventsSummary = repmat(1,m,3);
    for d = 1:size(unique(AllData(:,1)),1),      % first just do start/end times for events-not subdivided
        UEventIND = find(AllData(:,1) == d);    %index of all rows for event d, in AllData
        startE = AllData(UEventIND(1),4);
        endE = AllData(UEventIND(end),4);
        if d == 1,
            firstrow = [d startE endE];
            AllEventsSummary = firstrow;
        else
            nextrow = [d startE endE];
            AllEventsSummary = [AllEventsSummary; nextrow];
        end
    end
    
    % loop through AllEventsSummary, and remove any event that doesn't have
    % positive YDYN surge values reaching past a given threshold
    for s = 1:size(AllEventsSummary,1)
%         SdataIst = find(A(:,4) == AllEventsSummary(s,2));
%         SdataIend = find(A(:,4) == AllEventsSummary(s,2)+(timestep.*round(fs./2)));
        SdataI = find(A(:,4) >= AllEventsSummary(s,2) & A(:,4) <= AllEventsSummary(s,2)+(ratestep.*round(fs./2)));
        maxS = max(YDYN(SdataI));
        if maxS < surgeTH
            AllEventsSummary(s,1) = 0;
        end
    end
    
    ScutI = find(AllEventsSummary(:,1) == 0);
    AllEventsSummary(ScutI,:) = [];
    for ADloop = 1:length(ScutI)                            % remove these from AllData as well
        AllDataCUT = find(AllData(:,1)==ScutI(ADloop));
        AllData(AllDataCUT,1) = 0;
    end
    ADcutI = find(AllData(:,1) == 0);
    AllData(ADcutI,:) = [];
   
    
else disp('No Events Detected');
end
% so now column 7 has the within-event sub-events labelled. 

% %-----------------------------------------------------------------------------------------------------------------------
% % Plotting code (PART 1)
% 
% plotYthr = repmat(YDYNthreshold,1,size(AJdata,1));      % create YDYN threshold line
% plotZthr = repmat(ZJthreshold,1,size(AJdata,1));        % create Zjerk threshold line
% 
% %'time' as x-axis if present
% %timecheck = exist('time','var');
% if timecheck == 0
%     XVAR = time;
% elseif timecheck == 1
%     XVAR = AJdata(:,4); 
% end
% 
% %Plot the skeleton
% figure
% ax(1)=subplot(2,1,1);                                   % Subplot 1: YDYN
% h1c = plot(XVAR,YDYN);                           % plot YDYN vs time
% hold on
% h1d = plot(XVAR,plotYthr,'--');                  % plot YDYN threshold
% hold on
% ax(2)=subplot(2,1,2);                                   % Subplot 2: smoothzjdata vs time
% h2b = plot(XVAR,smoothzjdata);                   % smoothzjdata vs time
% hold on
% h2c = plot(XVAR,plotZthr,'--');                  % plot ZJ threshold
% hold on            
% linkaxes(ax,'x');
% 
% % plotting preferences
% set(h1c,'LineWidth',1.5);           % YDYN line width (subplot 1)
% set(h1c,'Color',[0.8 0.35 0.2]);    % YDYN line color
% set(h1d,'Color',[0.25 0.25 0.25]);  % YDYNthreshold color
% set(h2b,'Color',[0 0 0.3]);         % smoothzjdata color
% set(h2b,'LineWidth',1.5);           % smoothzjdata line width
% set(h2c,'Color',[0.4 0.4 0.4]);     % ZJthreshold line color
% 

%---------------------------------------------------------------------------------------------------------------------
% currently AllData(:,4) shows the 'created' timestamps --> if 'time' is an
% input, use an index to find the corresponding 'true' timestamp and use
% those instead.
if timecheck == 0
    realTime = ones(size(AllData,1),1);
    for ADlen = 1:size(AllData,1)
        tAD = AllData(ADlen,4);
        dif = abs(AJdata(:,4) - tAD);
        mD = min(dif);
        ADTind = find(dif == mD);
        if ~isempty(ADTind)
            ADTi = ADTind(1);
            realTime(ADlen) = time(ADTi);
        else
            realTime(ADlen) = 0;
        end
    end
%     AllData(:,4) = realTime;
end



% Apply depth threshold: if depth_threshold is added as an input arguement,
% remove events shallower than the threshold from AllData
if DThresh_check == 0
    Us = unique(AllData(:,1));
%     Dcut = ones(length(Us),1);
    for L_U = 1:size(Us,1)
        Iu = find(AllData(:,1) == Us(L_U));
        T1=Iu(1);
        Tu = realTime(T1);         % 1st accel timestamp for each event
        step = (time(2)-time(1))*(rateA/2);
        TrangeD = find(D_T > (Tu - step) & D_T < (Tu + step));
        meanD = mean(depth(TrangeD));
        if abs(meanD) < depth_threshold
            AllData(Iu,1) = 0;
        end
    end
    I_DTh = find(AllData(:,1) == 0);
    AllData(I_DTh,:) = [];
end
%-----------------------------------------------------------------------------------------------------------------------
% Plotting code (PART 1)

plotYthr = repmat(YDYNthreshold,1,size(AJdata,1));      % create YDYN threshold line
plotZthr = repmat(ZJthreshold,1,size(AJdata,1));        % create Zjerk threshold line

%'time' as x-axis if present
%timecheck = exist('time','var');
if timecheck == 0
    XVAR = time;
elseif timecheck == 1
    XVAR = AJdata(:,4); 
end

%Plot the skeleton
figure
ax(1)=subplot(2,1,1);                                   % Subplot 1: YDYN
h1c = plot(XVAR,YDYN);                           % plot YDYN vs time
hold on
h1d = plot(XVAR,plotYthr,'--');                  % plot YDYN threshold
hold on
ax(2)=subplot(2,1,2);                                   % Subplot 2: smoothzjdata vs time
h2b = plot(XVAR,smoothzjdata);                   % smoothzjdata vs time
hold on
h2c = plot(XVAR,plotZthr,'--');                  % plot ZJ threshold
hold on            
linkaxes(ax,'x');

% plotting preferences
set(h1c,'LineWidth',1.5);           % YDYN line width (subplot 1)
set(h1c,'Color',[0.8 0.35 0.2]);    % YDYN line color
set(h1d,'Color',[0.25 0.25 0.25]);  % YDYNthreshold color
set(h2b,'Color',[0 0 0.3]);         % smoothzjdata color
set(h2b,'LineWidth',1.5);           % smoothzjdata line width
set(h2c,'Color',[0.4 0.4 0.4]);     % ZJthreshold line color

%-------------------------------------------------------------------------------------------------------------------------
% create DATAevents: a data table quantifying and summarizing both full-event and sub-event (each motion within an event) data

if ~isempty(AllData),
    UniqueEvents = unique(AllData(:,1));
    DATAevents = ones(size(AllEventsSummary,1),7);
    
    for d = UniqueEvents(1):UniqueEvents(end),
        row = find(UniqueEvents == d);
        DATAevents(row,1) = row;                                                        % column 1: event number
        UEventIND = find(AllData(:,1) == d);
        if ~isempty(UEventIND)
            DATAevents(row,2) = AllData(UEventIND(end),4)-AllData(UEventIND(1),4);    % column (2 now): overall event duration
            DATAevents(row,3) = max(AllData(UEventIND(1:end),5));                     % column (3): max Z-axis Jerk
            DATAevents(row,4) = max(AllData(UEventIND(1:end),6));                     % column (4): max YDYN
            % add integral of ZJ and |YDYN|
            DATAevents(row,5) = trapz(AllData(UEventIND(1:end),5));                     % column 5: integral of ZJ
            DATAevents(row,6) = trapz(abs(AllData(UEventIND(1:end),6)));                % column 6: integral of absolute value of Ydyn
            %add start time to summary
            DATAevents(row,7) = AllData(UEventIND(1),4);

            %indicate start of feeding events with green squares
            plotData = AllData(UEventIND(1),4);
            if timecheck == 0
                PDind = find(A(:,4) == plotData);     %x axis is 'time' variable, if given as input
                Xtime = XVAR(PDind);
            else
                Xtime = plotData;
            end
            plotX = repmat(0,size(plotData,1),1);
            subplot(2,1,1)
            plot(Xtime,plotX,'gs')
            hold on
        end
    end
    summary = DATAevents;
    S_cut = find(summary(:,1) == 1 & summary(:,2) == 1 & summary(:,3) == 1);
    summary(S_cut,:) = [];
    Sum_header = {'Event_#', 'Duration', 'Max_Heave_Jerk', 'Max_Surge_Decel', 'Heave_Jerk_Int', 'Abs_Surge_Decel_Int', 'Start_Time'};
    Summary = [Sum_header; num2cell(summary)];
else disp('No Feeding Detected')    
end

Event_Data = [AllData(:,1) AllData(:,4:6)];

listE = unique(Event_Data(:,1));
for L_E = 1:size(listE,1)
    Elist_E = find(Event_Data(:,1) == listE(L_E));
    Event_Data(Elist_E,1) = L_E;
end

%summary = DATAevents;
% Sum_header = {'Event_#', 'Duration', 'Max_Heave_Jerk', 'Max_Surge_Decel', 'Heave_Jerk_Int', 'Abs_Surge_Decel_Int', 'Start_Time'};
% Summary = [Sum_header; num2cell(summary)];

%2nd figure: depth with feeding events, ONLY if depth and rateD are both added as inputs
depth_check = isempty(depth);
if depth_check == 0
    DT_check = isempty(depth_time);
    if DT_check == 0
        DT = depth_time;
        Tstep_d = DT(2) - DT(1);
        tDetec = ones(size(summary,1),1);
        for dets = 1:size(summary,1)
            %first, find where A(:,4) == summary(:,7) -> use as index of where
            %events are in depth_time
            IDet = find(A(:,4) == summary(dets,7));
            IDet = IDet(1);
            det_DT_est = round(IDet./(rateA/rateD));               %estimate of where detections should be in the depth_time vector
            tDetec(dets) = DT(det_DT_est);
            rDet = find(DT > (tDetec(dets) - (Tstep_d.*10)) & DT < (tDetec(dets) + (Tstep_d.*10)));
            DetDep = depth(rDet);
            uDep = mean(DetDep);
            dDetec(dets) = uDep;
        end
    elseif DT_check == 1
        if timecheck == 1
            %start depth_time vector at 0 with a step length of 'rateD'
            rateDcheck = exist('rateD','var');
            if rateDcheck == 0
                DT = (0:(1/rateD):XVAR(end));
                DT = DT(1:length(depth));
                tDetec = summary(:,7)./(rateA/rateD);
                for L_dep = 1:size(summary,1)
                    rDet = find(DT > (tDetec - (rateD.*10)) & DT < (tDetec + (rateD.*10)));
                    DetDep = depth(rDet);
                    uDep = mean(DetDep);
                    dDetec(dets) = uDep;
                end
            else
                disp('input rateD to plot detected feeding with depth')
            end
        elseif timecheck == 0
            %start depth_time vector at time(1), step length of 'timestep'
            timestep = time(2)-time(1);
            rateDcheck = exist('rateD','var');
            if rateDcheck == 0
                rateFac = (rateA/rateD);
                DT = (time(1):(timestep/rateFac):XVAR(end));            %estimated depth_time vector using inputed rates
                DT = DT(1:length(depth));
                for dets = 1:size(summary,1)
                    %first, find where A(:,4) == summary(:,7) -> use as index of where
                    %events are in depth_time
                    IDet = find(A(:,4) == summary(dets,7));
                    det_DT_est = round(IDet./(rateA/rateD));               %estimate of where detections should be in the depth_time vector
                    tDetec(dets) = DT(det_DT_est);
                    rDet = find(DT > (tDetec(dets) - ((timestep/rateFac).*10)) & DT < (tDetec(dets) + ((timestep/rateFac).*10)));
                    DetDep = depth(rDet);
                    uDep = mean(DetDep);
                    dDetec(dets) = uDep;
                end
            else
                disp('input rateD to plot detected feeding with depth')
            end
        end
    end
end

if depth_check == 0
    figure
    plot(DT,depth,'LineWidth',1,'Color',[0 0 0]);
    hold on
    plot(tDetec,dDetec,'ro');
    hold on
end

end
%-------------------------------------------------------------------------------------------------------------------------------
