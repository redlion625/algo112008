function gast = cal2GAST(y,m,d,hh,mm,ss)
% Given:
%     y,md,hh,mm,ss in UTC 
%
% Returned:
%     gast    Greenwich apparent sidereal time in radians
% 
utc=greg2mjd(y,m,d,hh,mm,ss);
uta = Const.DJM0;
[DUT1] = extractDUT1(y,m,d);
utb=utc+DUT1/Const.DAYSEC;

tta=Const.DJM0;
ttb=utc+(Const.TTMTAI + Const.TAIUTC)/Const.DAYSEC;

gast = gst06a(uta,utb,tta,ttb);





