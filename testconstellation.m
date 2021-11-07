function output = globalGNSS(line)
str = 'G C E R';
expression = line;
startIndex = regexp(str,expression);
if(startIndex > 0)
    fprintf("This is correct!\n"); 
else
    fprintf("This is wrong!\n"); 
end

line = 'C'; 

% str = 'bat cat can car coat court CUT ct CAT-scan';
% expression = 'c[aeiou]+t';
% startIndex = regexp(str,expression)

str = 'G C E R';
expression = line;
startIndex = regexp(str,expression);
if(startIndex > 0)
    fprintf("This is correct!\n"); 
else
    fprintf("This is wrong!\n"); 
end
