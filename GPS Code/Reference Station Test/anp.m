function w = anp(a)
w = rem(a, Const.D2PI);

if (w < 0)
    w = w + Const.D2PI;
end