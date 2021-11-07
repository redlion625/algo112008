function str = fixedWidth(line,widths)
% fixedWidth(line,width)
% Takes string or vector char array and integer vector of widths to separate
% said line into string array. Widths are the width of the 
% separated field starting and accumulating from the left of the line
%
% sum(widths) <= numel(line)
%
% For example
% 
% line='G06  22327780.74706   1190616.93306     -2779.59406        36.90006
% 22325439.968 3    877396.458 3     -2074.980 3        20.700 3';
% 
% widths=[3 16 16 16 16 14 2 14 2 14 2 14 2];
%
% str=fixedWidth(line,widths)
% 
% str = 
% 
%   1×13 string array
% 
%   Columns 1 through 9
% 
%     "G06"    "  22327780.74706"    "   1190616.93306"    "     -2779.59406"    
%     "        36.90006"    "  22325439.968"    " 3"    "    877396.458"    " 3"
% 
%   Columns 10 through 13
% 
%     "     -2074.980"    " 3"    "        20.700"    " 3"
% See also str2num
N=numel(widths);
str = strings(1,N);
if isstring(line)
    line=char(line);
end

if sum(widths)>numel(line)
    error('Sum of widths exceed length of line, empty string array returned')
else
    widths=reshape(widths,1,N);
    widths=[0 widths];
    for i=1:N
        startInd=sum(widths(1:i))+1;
        endInd=sum(widths(1:i+1));
        str(i)=line(startInd:endInd);
    end
end
end