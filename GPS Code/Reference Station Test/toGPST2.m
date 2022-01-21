function gps_sec=toGPST2(y,m,d,hh,mm,ss)
% Converts Gregorian Calendar date in GPS time (no leap seconds) to number of GPS seconds
% since GPS start epoch Jan. 6, 1980 00:00:00
gps_sec=secondsSince(y,m,d,hh,mm,ss,1980,1,6,0,0,0,0);
%gps_sec=mod(gps_sec,604800)% seconds of week, assumes week starts on Sunday 00:00:00