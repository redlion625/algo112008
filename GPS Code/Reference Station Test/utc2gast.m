function gast = utc2gast(y,m,d,hh,mm,ss)
% utc2gast(y,m,d,hh,mm,ss)
% Given:
%     y,m,d,hh,mm,ss in UTC 
%
% Returned:
%    Greenwich apparent sidereal time in radians
utc=greg2mjd(y,m,d,hh,mm,ss);
uta = Const.DJM0;
% DUT1 = deltaUT1(utc);
addpath(genpath('lib'));
[DUT1,DUT1err] = extractDUT1(y,m,d);
utb=utc+DUT1/Const.DAYSEC;

tta=Const.DJM0;
ttb=utc+(Const.TAIUTC+Const.TTMTAI)/Const.DAYSEC;

gast = gst06a(uta,utb,tta,ttb);





