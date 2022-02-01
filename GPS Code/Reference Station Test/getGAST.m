function GAST = getGAST(jd,DUT1,DTAI)%,outType)
% getGAST - uses ERA(UT1) based computation of Greenwich Apparent Sidereal
% Time (GAST) based on IERS 2010 Conventions
%Inputs:
%jd - Julian date in UTC at GMT
%DUT1 - (UT1-UTC) from IERS Bulletin A in seconds
%DTAI - (TAI-UTC) in seconds
% e.g. for 04/28/2021
%DTAI=37s from Jan. 2017 - ???
%DUT1=-0.182695 for 04/28/2021 from IERS Bulletin A
%
%Output:
% GAST in hour angle (DMS)

jdTT=jd+(DTAI+32.184)/86400; % UTC to TAI to TT
t=(jdTT-2451545)/36525; % TT Julian date since J2000

T_u=(jd+DUT1)-2451545; % UT1 Julian date since J2000

EOpoly=0.014506 + 4612.156534*t + 1.3915817*t^2 ...
    - 0.00000044*t^3 - 0.000029956*t^4 - 0.0000000368*t^5; % in arcsecond



% Earth
% ERA=2*pi*(0.7790572732640 + 1.00273781191135448*T_u);

% equivalent ERA calculation robust to rounding errors
jdUT1frac=T_u-floor(T_u); % UT1 fraction of day
ERA=2*pi*(jdUT1frac+0.7790572732640 + 0.00273781191135448*T_u);

% T=(jd-2451545)/36525;
% GMST = 24110.54841 + 8640184.812866*T + 0.093104*T^2 - 6.2e-6*T^3 % in seconds of hour angle?
%
% GMST=mod(GMST,86400); % wrap to [0h,24h]
% if strcmp(outType,'hms')
%     hh = floor(GAST/3600);
%     mm = floor((GAST-3600*hh)/60);
%     ss = GAST-3600*hh-60*mm;
%     GAST = [hh mm ss];
% elseif strcmp(outType,'rad')
%     GAST=2*pi*GAST/86400;
% elseif strcmp(outType,'deg')
%     GAST=GAST/240;
% else
%     error('outType not recognized')
% end
