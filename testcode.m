function str = spacereplace(line)
%line = '33'; 
%line = '28497045.494 3';
    testarray = find((isspace(line)==1)); 
    inttest = int64(testarray);
    if(inttest > 0)
        str = line; 
        expression = ' '; 
        replace = '0'; 
        line = regexprep(str,expression,replace);
    else
        line = line; 
    end
    again = 3; 
end
str = line; 

% line(2:3) = '33'; 
% line(70:83) = '28497045.494 3'; 
% 
% test = line(70:83);  

% testarray = find((isspace(line(70:83))==1)); 
% %a,'uint8')
% inttest = int64(testarray); %cast(testarray,'uint8');% integer(testarray);
% if(isinteger(inttest) ==1)
% %if(isdouble(testarray) == 1)
%     str = line(70:83); 
%     expression = ' '; 
%     replace = '0'; 
%     line(70:83) = regexprep(str,expression,replace);
% %     test(inttest) = char(1);
% %     line(70:83) = test; 
% end
% % array = isspace(line(2:3)); 
% % if((find(array)==1) == 
% 

