function gast = cal2GAST(y,m,d,hh,mm,ss,DUT1)
% Given:
%     y,md,hh,mm,ss in UTC 
%     DUT1 of date (UT1-UTC)
% 
utc=greg2mjd(y,m,d,hh,mm,ss);
uta = Const.DJM0;
utb=UTC+DUT1/Const.DAYSEC;

tta=Const.DJM0;
ttb=UTC+(Const.TAIUTC+Const.TTMTAI)/Const.DAYSEC;

gast = gst06a(uta,utb,tta,ttb);





