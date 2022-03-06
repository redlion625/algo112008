%% toGPST
% The following function converts the given epoch to BDS
% Input:    year,month,day,hour,minute,second
% Output:   BDS week and BDS seconds of week

function gps_seconds = toBDS(year,month,day,hour,minute,second)
    if month <= 2
        year = year - 1;
        month = month +12;
    end
    hour = hour + minute/60 + second/3600; %real hour
    jd = floor(365.25*year) + floor(30.6001*(month+1)) + day + hour/24 + 1720981.5; %convert to Julian Date
    week = floor((jd - 2453736.5)/7); %compute BDS week 
    frc_of_day = mod(floor(jd+0.5),7) + jd - floor(jd) - 0.5; % Day of week + fraction of day - noon
    if hour >= 12
        frc_of_day = frc_of_day+1;
    end
    secondsofweek = (frc_of_day+1)*86400; % Add one since BDS weeks start at zero
    gps_seconds = secondsofweek + week*60*60*24*7; %BDS seconds since 1-Jan-2006
end