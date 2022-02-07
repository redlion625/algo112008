function jd = greg2jd(y,m,d,hh,mm,ss)
% greg2jd - time in decimal days since noon at Greenwich at -4713 November
% 23
% e.g. J2000 epoch greg2jd(2000,1,1,12,0,0) = 2451545

% Fliegel and Van Flandern algorithm for converting
% Gregorian calendar date to julian day number jd

%SoD=hh*3600+mm*60+ss
%jd=(1461*(y+4800+(m-14)/12))/4+(367*(m-2-12*((m-14)/12)))/12 ...
%    -(3*((y+4900+(m-14)/12)/100))/4+d-32075;%+SoD/86400;

% full form algorithm *works sometimes
%jd=367*y-7*floor((y+floor((m+9)/12))/4) ... 
%    -3*floor(floor((y+floor((m-9)/7))/100)/4) ...
%    +floor(275*m/9)+d+1721029;

% short form algorithm only valid for dates after March 1900
jd=367*y-floor(7*(y+floor((m+9)/12))/4)+floor(275*m/9)+d+1721014;

jd=jd+(hh*3600+mm*60+ss-43200)/86400; % fraction of day