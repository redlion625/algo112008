function [GMST] = getGMST(jd,outType)
% getGMST - calculates GMST, eastward angle between vernal equinox and
%           and Greenwich Meridean ignoring variations from nutation 
% Input:    jd - julian date
%           outType - string or char array of output type
%                'hms' for hh mm ss hour angle
%                'rad' for radians
%                'deg' for decimal degrees
% Output: Greenwich Mean Sidereal Time in radians

T=(jd-2451545)/36525;
GMST = 24110.54841 + 8640184.812866*T + 0.093104*T^2 - 6.2e-6*T^3 % in seconds of hour angle?
GMST=mod(GMST,86400); % wrap to [0h,24h]
if strcmp(outType,'hms')
    hh = floor(GMST/3600);
    mm = floor((GMST-3600*hh)/60);
    ss = GMST-3600*hh-60*mm;
    GMST = [hh mm ss];
elseif strcmp(outType,'rad')
    GMST=2*pi*GMST/86400;
elseif strcmp(outType,'deg')
    GMST=GMST/240;
else
    error('outType not recognized')
end
    