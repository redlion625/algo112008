function [sod,sow,gpsweek] = toGPST3(y,m,d,hh,mm,ss)
gpsutc=18;
%y=2021;m=4;d=28;hh=22;mm=19;ss=8;%hh=0;mm=0;ss=0; 
% expected result
%[ 80366,339566,2155]
jd=greg2jd(y,m,d,0,0,0);
jdgps=jd-2444244.5;
gpsweek=floor(jdgps/7);
dow=jdgps-gpsweek*7; %day of week from 0-6 starting at 0 = Sunday;
sod=hh*3600+mm*60+ss+gpsutc;
sow=dow*Const.DAYSEC+sod;


