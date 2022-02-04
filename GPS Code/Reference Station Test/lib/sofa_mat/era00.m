function theta = era00(dj1,dj2)
if (dj1 < dj2)
    d1 = dj1;
    d2 = dj2;
else
    d1 = dj2;
    d2 = dj1;
end
t = d1 + (d2 - Const.DJ00);

% Fractional part of T (days).
f = rem(d1, 1.0) + rem(d2, 1.0);

% Earth rotation angle at this UT1.
theta = anp(Const.D2PI * (f + 0.7790572732640 ...
+ 0.00273781191135448 * t));