function gast = gst06(uta,utb,tta,ttb,rnpb)

[x, y] = bpn2xy(rnpb);

% The CIO locator, s. */
s = s06(tta, ttb, x, y);

% Greenwich apparent sidereal time. */
era = era00(uta, utb);
eo = eors(rnpb, s);
gast = anp(era - eo);