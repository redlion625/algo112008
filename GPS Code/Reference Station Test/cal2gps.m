function t = cal2gps(y,m,d,hh,mm,ss)
% cal2gps - converts calendar date to GPS time in seconds
%
% Given:
%     y,m,d,hh,mm,ss in gps time
%     
% Returned:
%     t     number of seconds since GPS start epoch January 6, 1980 0h
%


jd = greg2jd(y,m,d,0,0,0); 
jdgps = jd-2444244.5; 
t = jdgps*86400 + 3600*hh + 60*mm + ss;
