function t = utc2gps(y,m,d,hh,mm,ss)
% utc2gps - converts calendar date in utc to gps time
%
% Given:
%     y,m,d,hh,mm,ss in utc
%     
% Returned:
%     t     number of seconds since GPS start epoch January 6, 1980 0h
%
% Notes:
%         Only applies to dates from 2017/01/01 - 2022/12/31 due to GPS-UTC
%         leap second correction assumption

gpsutc = 18; %(GPST - UTC)

% jd = greg2jd(y,m,d,0,0,0); 
% jdgps = jd-2444244.5; 
% t = jdgps*86400 + 3600*hh + 60*mm + ss + gpsutc;

t = cal2gps(y,m,d,hh,mm,ss) + gpsutc;