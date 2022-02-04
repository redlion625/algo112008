function era = eraGLO(jd)
% eraGLO(jd)
era = 2*pi*(0.7790572732640 + 1.00273781191135448*(jd-Const.DJ00));