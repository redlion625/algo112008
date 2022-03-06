function obs = obs_read_rinex211(filename)
tic
fclose('all');
%filename = 'algo1180.21o'
fid=fopen(filename);
line=fgetl(fid);
obs=struct();
obs.obstypes = [];
obs.interval = [];
while ~contains(line,'END OF HEADER')
    if contains(line,'TYPES OF OBS')
        linesplit = fixedWidth(line,6*ones(1,9));
        obs.obstypes = strip(linesplit(linesplit~=""));
    elseif contains(line,'Interval')
        obs.interval = num2str(strip(line(1:60)));
    elseif contains(line,'TIME OF FIRST OBS')
        obs.date_obs = str2double(fixedWidth(line,[6 6 6]));
    end
    line = fgetl(fid);
end
if isempty(obs.interval)
    epoch = cell(1,3);
else
    epoch = cell(Const.DAYSEC/obs.interval,2);
end

line = fgetl(fid);
pat = lettersPattern;
count = 0;
% while end of file not found
while ~feof(fid)
    if contains(line,pat)
        count = count + 1;
        linesplit = fixedWidth(line,[3*ones(1,5) 11 3 3 36]);
        ymdhms = str2double(linesplit(1:6));
        numsat = str2num(linesplit(8));
        
        epoch{count,3}=ymdhms; % calendar time GPST
        epoch{count,4}=ymdhms(4)*3600+ymdhms(5)*60+ymdhms(6); % second of day GPST
        
        epoch{count,1} = cal2gps(2000+ymdhms(1),ymdhms(2),...
            ymdhms(3),ymdhms(4),ymdhms(5),ymdhms(6));
        obsrec = NaN(numsat,10);
        linesplit = linesplit(9);
        if numsat > 12
            line = fgetl(fid);
            linesplit = strcat(linesplit,strip(line));
        end
        widths = zeros(1,numsat*2);
        widths(1:2:end) = 1;
        widths(2:2:end) = 2;
        linesplit = fixedWidth(linesplit,widths);
        sys = linesplit(1:2:end);
        
        % fill 2nd col with PRN
        obsrec(:,2) = str2double(linesplit(2:2:end)');
        for i = 1:numsat
            % fill 1st column with constellation codes
            % G = 1, R = 2, else = 69
            obsrec(i,1) = const2num(sys(i));
            
            % read in dual frequency observables with headers
            % stored in obs.obstypes
            line = fgetl(fid);
            linesplit = fixedWidth(line,[14 16 17 16 15]);
            
            line = fgetl(fid);
            linesplit = cat(2,linesplit,fixedWidth(line,[14 16 17]));
            
            obsrec(i,3:10) = str2double(linesplit);
        end
        
        epoch{count,2} = obsrec;
                
        line = fgetl(fid);
    else
        line = fgetl(fid);
        error('what is this file?');
    end
    
    obs.epochdata = epoch;
end
toc



