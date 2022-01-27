
function secSince = secondsSince(y2,m2,d2,hh2,mm2,ss2,y1,m1,d1,hh1,mm1,ss1,leapsec)
% secondsSince counts seconds since epoch not including leap seconds given Gregorian
% calendar date
% Used to calculate seconds since start epoch1 of any continuous or atomic time scales
% given the right start epoch until end epoch2
% Notes:    does not consider day added every 400 years

dpm=[0 31 28 31 30 31 30 31 31 30 31 30 31]; % days per month

dayFrac1=(60-ss1)+(59-mm1)*60+(23-hh1)*3600; % in seconds
dayFrac2=ss2+mm2*60+3600*hh2;

% if day of year before or after midnight before potential leap day

doy1=sum(dpm(1:m1))+d1;

if doy1<60
    dToFeb28=60-doy1-1;
    if mod(y1,4)==0
        y2leap=0;
    else
        y2leap=4-mod(y1,4);
    end
    leap1=y1+y2leap;
else
    dToFeb28=435-doy1;
    if mod(y1+1,4)==0
        y2leap=0;
    else
        y2leap=4-mod(y1+1,4);
    end
    leap1=y1+y2leap+1;
end
dToLeap=dToFeb28+y2leap*365;
doy2=sum(dpm(1:m2))+d2;

if doy2<60
    feb28ToD=305+doy2;
    leap2y=mod(y2-1,4);
    leap2=y2-leap2y-1;
else
    feb28ToD=doy2-60;
    leap2y=mod(y2,4);
    leap2=y2-leap2y;
end
numLeapDays=(leap2-leap1)/4+1;

leapToD=feb28ToD+leap2y*365;

secSince=dayFrac1+dToLeap*86400+(leap2-leap1)*31536000+numLeapDays*86400+leapToD*86400+dayFrac2+leapsec;
end
