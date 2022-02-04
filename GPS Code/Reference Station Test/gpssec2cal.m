function cal = gpssec2cal(gpssec)
% Given: GPS seconds since GPS start epoch January 6, 1980 00:00:00
% 
% Returned: date vector in format [year month day hour min sec] in GPST
%

days = floor(gpssec/Const.DAYSEC);

jd = 2444244.5 + days;
[y,m,d] = jd2greg(jd);

sod = gpssec - days*Const.DAYSEC;
hh = floor(sod/3600);
mm = floor((sod - hh*3600)/60);
ss = sod - hh*3600 - mm*60;

cal = [y m d hh mm ss];


