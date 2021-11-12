function nav = navreadGLONASS(filename)
% Reads navigation file from filename
% from GLONASS of RINEX v. 2 format
% into nav struct with an instance for each epoch

%filename='brdc1180.21g';
fid=fopen(filename);

% skip irrelevant info lines
for i = 1:4
    fgetl(fid);
end

% store useful header info
line=fgetl(fid);
temp=fixedWidth(line,[6 6 6 22]);
year=str2num(temp(1));
month=str2num(temp(2));
day=str2num(temp(3));
glonass2utc=str2num(temp(4));

line=fgetl(fid);
leap=str2num(fixedWidth(line,[6]));
fgetl(fid);
nav=struct([]);

% end of header reached navigation mess begin
i=1;
line = fgetl(fid);
while ~feof(fid)

    % give nav msg struct per epoch i header info
    nav(i).year=year;
    nav(i).month=month;
    nav(i).day=day;
    nav(i).glonass2utc=glonass2utc;

    % read line 1 of mess at epoch i
    temp=fixedWidth(line,[2 9 3 3 5 19 19 19]);
    nav(i).satnum=str2num(temp(1));
    nav(i).hr=str2num(temp(3));
    nav(i).min=str2num(temp(4));
    nav(i).sec=str2num(temp(5));
    nav(i).clkbias=str2num(temp(6));
    nav(i).freqbias=str2num(temp(7));
    nav(i).msgframetime=str2num(temp(8));

    % read line 2
    line=fgetl(fid);
    temp=fixedWidth(line,[22 19 19 19]);
    nav(i).xpos=str2num(temp(1));
    nav(i).xvel=str2num(temp(2));
    nav(i).xacc=str2num(temp(3));
    nav(i).sathealth=str2num(temp(4));

    % read line 3
    line=fgetl(fid);
    temp=fixedWidth(line,[22 19 19 19]);
    nav(i).ypos=str2num(temp(1));
    nav(i).yvel=str2num(temp(2));
    nav(i).yacc=str2num(temp(3));
    nav(i).freqnum=str2num(temp(4));

    % read line 4
    line=fgetl(fid);
    temp=fixedWidth(line,[22 19 19 19]);
    nav(i).zpos=str2num(temp(1));
    nav(i).zvel=str2num(temp(2));
    nav(i).zacc=str2num(temp(3));
    nav(i).ageop=str2num(temp(4));

    i=i+1;
    line=fgetl(fid);

    % for testing
    %     if i>1
    %         break
    %     end
end

end
% line = fgetl(fid);
% end_found = contains(line,'END OF HEADER'); % search to find end of header

% while ~end_found %skips header lines
%     line=fgetl(fidobs);
%     end_found = contains(line,'END OF HEADER'); % search to find end of header
%
% end