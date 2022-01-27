function GAST = getGAST(jd,outType)
T=(jd-2451545)/36525;
GMST = 24110.54841 + 8640184.812866*T + 0.093104*T^2 - 6.2e-6*T^3 % in seconds of hour angle?

GMST=mod(GMST,86400); % wrap to [0h,24h]
if strcmp(outType,'hms')
    hh = floor(GAST/3600);
    mm = floor((GAST-3600*hh)/60);
    ss = GAST-3600*hh-60*mm;
    GAST = [hh mm ss];
elseif strcmp(outType,'rad')
    GAST=2*pi*GAST/86400;
elseif strcmp(outType,'deg')
    GAST=GAST/240;
else
    error('outType not recognized')
end
