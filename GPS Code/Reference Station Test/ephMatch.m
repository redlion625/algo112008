function minIndex = ephMatch(satnumobs,tobs,satnumnav,tnav)
% ephMatch(satnumobs,tobs,satnumnav,tnav)
% Given:
%     1 x 1 double        satnumobs       PRN or satellite number of
%                                         observation 
%     1 x 1 double        tobs            time of observation 
%     n x 1 double        satnumnav       satellite numbers of a list of n
%                                         navigation records
%     n x 1 double        tnav            time of navigation records
%     
% Returned:
%     1 x 1 double        minIndex        index of navigation record
%                                         leading to the minimum time 
%                                           difference
%                                     

indexMatch = [];
timeMatch = [];

%% Filters for matching mavigation records by PRN

% Searches for PRNs in each navigation record matching the observation PRN
for i = 1:numel(satnumnav)
    
    % store navigation PRNs matching the observation PRN and the navigation
    % record's index
    if satnumobs == satnumnav(i)
        
        indexMatch = [indexMatch;i];
        timeMatch = [timeMatch, tnav(i)];
    end
end

%% Find minimum time difference

% Initialize minimum time difference
minDeltaT = abs(tobs-timeMatch(1));
minIndex = indexMatch(1);

for i = 2:numel(indexMatch)
    
    deltaT = abs(tobs-timeMatch(i));
    
    % checks each navigation record's time, replaces the minimum time
    % difference whenever a new minimum is detected and stores the
    % corresponding navigation record's index
    if deltaT < minDeltaT
        
        minDeltaT = deltaT;
        minIndex = indexMatch(i);
    end
end
    
        