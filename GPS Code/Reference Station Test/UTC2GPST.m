function gps_sec=UTC2GPST(y,m,d,hh,mm,ss,leapsec)
% Converts Gregorian Calendar date in UTC to number of GPS seconds
% since GPS start epoch Jan. 6, 1980 00:00:00 given correct number of leap
% seconds (GPS-UTC) given
gps_sec=secondsSince(y,m,d,hh,mm,ss,1980,1,6,0,0,0,leapsec);