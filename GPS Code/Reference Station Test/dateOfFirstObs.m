function date = dateOfFirstObs(filename)
% Searches RINEX 3.04 observation file for date of first observation
%
% Given: filename of observation file
% 
% Returned: date vector as [ year month date ]


fid = fopen(filename);

line = fgetl(fid);
while ~contains(line, 'END OF HEADER')
    if contains(line, 'TIME OF FIRST OBS')
        
        % First 3 entries are year, month, day
        date = str2double(fixedWidth(line,[6 6 6]));
        return
        
    else
        line = fgetl(fid);
    end
end
warning('Date not found, -1 returned');
date = -1;