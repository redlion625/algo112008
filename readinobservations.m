%% Observation File Reader
%Updated notes: The first column of the observation file is the
%constellation (GPS = 1, GLONASS = 2, GALILEO = 3), the second column is
%the satellite number, the third column is the pseudorange in Frequency 1,
%the fourth column is the phase in Frequency 1, the fifth column is the
%doppler in Frequency 1, the sixth column is the raw signal strength in
%Freequency 1, the seventh column is the pseudorange in Frequency 5, the
%eighth column is the phase in Frequency 5, the ninth column is the doppler
%in Frequency 5, and the tenth column is the raw signal strength in
%Frequency 5. 

% The following function extracts observation data from .##o files
% Input:    Observation filepath
% Output:   obs.rec_ID - Receiver ID
%           obs.pos_xyz - Estimated receiver position
%           obs.ant_delta_HEN - antenna position offset
%           obs.num_typeobs- number of observables (L1,L2,C1...)
%           obs.type_obs- list of observables (L1,L2,C1...)
%           obs.data_headers- list of observables stored in the data structure(L1,L2,C1...)
%           obs.data - corresponding data values for the variables stored
%           in the headers. Each cell row contains epoch information i.e. 
%           GPST and all obserables from satellites visible at the time.
%           Col 1 store time information, Col 2 stores satellite
%           information; PRN-L1-C1 so on. The observables are ordered per
%           the obs file layout

function obs = obs_read(filepath)
%D:\Third Year\ESSE 3670\Project 3\Data Downloaded\datasets to try with\algo\algo112008\ALGO0010.08O
file = fopen("Pixel4_GnssLog.21o"); %("Pixel4XLModded_GnssLog.20o");  %fopen("ALGO0010.08O");
disp('------------------Begin reading obs file---------------------');
%read constant parameters
curr_line = fgetl(file); %fget1 %fopen %fgetl
obs.rec_ID = []; 
obs.pos_xyz = [];
obs.ant_delta_HEN = [];
obs.num_typeobs = [];
obs.type_obs = [];


while isempty(strfind(curr_line,'END OF HEADER'))
    if ~isempty(strfind(curr_line,'REC # / TYPE / VERS'))
        temp = strsplit(curr_line);
        obs.rec_ID = temp(1);
    elseif ~isempty(strfind(curr_line,'APPROX POSITION XYZ'))
        temp = strsplit(curr_line);
        for j=1:3
            obs.pos_xyz(j) = str2double([temp{j+1}]);
        end
    elseif ~isempty(strfind(curr_line,'ANTENNA: DELTA H/E/N'))
        temp = strsplit(curr_line);
        for j=1:3
            obs.ant_delta_HEN(j) = str2double([temp{j+1}]);
        end
    elseif ~isempty(strfind(curr_line,'# / TYPES OF OBSERV'))
        temp = strsplit(curr_line);
        obs.num_typeobs = str2double(temp(2));
        obs.type_obs = temp(3:obs.num_typeobs+2);
    end
    curr_line = fgetl(file); %fgetl or fget1
end

%% Helper function
%  Standardizes length of data lines to dicern which observables are missing
function lineOut = fillWhite(line)
    lineOut = line;
   while length(lineOut) < 80
       
       lineOut = [lineOut,' '];
   end
    

end

obs.data_headers = {'Pseudorange', 'Phase', 'Doppler', 'Raw Signal Strength'};  %{'PRN',obs.type_obs{:},'Nav Ind','Sat ECF X','Sat ECF Y','Sat ECF Z','Delta tsv'}; 
obs.struct_headers = {'Pseudorange', 'Phase', 'Doppler', 'Raw Signal Strength'}; %{'Epoch Time','Data','Rec. ECF X','Rec. ECF Y','Rec. ECF Z','Rec. Clock Error','LS Iterations','Residuals','Cov. Matrix'};
data = cell(2880,2); % using 2 obs/min * 60 mins/hr 24hr/day

epoch = 1;
n = 1;
curr_line = fgetl(file);
while ~feof(file)
    % iterate per epoch
    % first line of data
    v = n;
    temp = strsplit(curr_line); 
    year = str2num(temp{2})-2000; %Epoch Year
    month = str2num(temp{3}); %Epoch Month
    day = str2num(temp{4}); %Epoch Day
    hour = str2num(temp{5}); %Epoch Hour
    minute = str2num(temp{6}); %Epoch Minute
    second = str2num(temp{7}); %Epoch Second
    flag = str2num(temp{8}); %Epoch Flag
    num_sat = str2num(curr_line(34:35)); %Number of measurements in epoch
    curr_line = fgetl(file);
    sat_names = curr_line(:,1:3); % Sat names
    sat_names= sat_names(find(~isspace(sat_names))); %remove empty spaces
    %create array of sat PRN's without 'G'
    sat_names = strsplit(sat_names,'G');
    sat_names = sat_names(~cellfun('isempty',sat_names));
    sat_PRN = zeros(num_sat,1);
    %store epoch time. CONVERT TO GPST and more
    data{epoch,1} = toGPST(year+2000,month,day,hour,minute,second);
    % store observables for each satellite 
    epoch_data = zeros(num_sat,obs.num_typeobs+6);
    epoch_data(:,1) = sat_PRN;
    line = curr_line; 
    for j = 1:num_sat
            str = 'G C E R';
            expression = line(1);
            startIndex = regexp(str,expression);
            if(startIndex > 0)
                if(line(1) == 'G')
                    epoch_data(j,1) = 1; 
                elseif(line(1) == 'R')
                    epoch_data(j,1) = 2; 
                elseif(line(1) == 'E')
                    epoch_data(j,1) = 3; 
                else
                    epoch_data(j,1) = 4; 
                end
                epoch_data(j,2) = str2num(line(2:3));
                
                if(83 > (numel(line))) 
                      epoch_data(j,7) = 0.0000000;            
                elseif((line(70:83)) == '              ')  
                      epoch_data(j,7) = 0.0000000;
                else
                      line(70:83) = spacereplace(line(70:83)); 
                      epoch_data(j,7) = str2num(line(70:83)); 
                end
                if(99 > (numel(line)))
                      epoch_data(j,8) = 0.0000000; 
                elseif((line(88:99)) == '            ')   
                      epoch_data(j,8) = 0.0000000;
                else
                      line(88:99) = spacereplace(line(88:99)); 
                      epoch_data(j,8) = str2num(line(88:99));
                end
                if(115 > (numel(line)))
                      epoch_data(j,9) = 0.0000000; 
                elseif((line(105:115)) == '           ')    
                      epoch_data(j,9) = 0.0000000;
                else
                     line(105:115) = spacereplace(line(105:115));
                     epoch_data(j,9) = str2num(line(105:115)); 
                end
                if(131 > (numel(line)))
                     epoch_data(j,10) = 0.0000000;
                elseif((line(124:131)) == '        ') 
                     epoch_data(j,10) = 0.0000000;
                else
                     line(124:131) = spacereplace(line(124:131));
                     epoch_data(j,10) = str2num(line(124:131));
                end
                
                %end
                if(19 > (numel(line)))
                   epoch_data(j,3) = 0.0000000; 
                elseif((line(6:19)) == '              ') 
                   %continue;
                   epoch_data(j,:) = 0.0000000;
                   line = fgetl(file);
                   continue;
                else
                   line(6:19) = spacereplace(line(6:19)); 
                   epoch_data(j,3) = str2num(line(6:19)); 
                end
                if(35 > (numel(line)))
                    epoch_data(j,4) = 0.0000000; 
                elseif((line(23:35)) == '             ')  
                    epoch_data(j,:) = 0.0000000;
                    line = fgetl(file);
                    continue;
                else
                     line(23:35) = spacereplace(line(23:35)); 
                     epoch_data(j,4) = str2num(line(23:35)); 
                end
                if(51 > (numel(line)))
                    epoch_data(j,5) = 0.0000000; 
                elseif((line(41:51)) == '           ')   
                    epoch_data(j,:) = 0.0000000;
                    line = fgetl(file);
                    continue; 
                else
                     line(41:51) = spacereplace(line(41:51)); 
                     epoch_data(j,5) = str2num(line(41:51)); 
                end
                if(67 > (numel(line)))
                     epoch_data(j,6) = 0.0000000; 
                elseif((line(60:67)) == '        ')  
                     epoch_data(j,:) = 0.0000000;
                     line = fgetl(file);
                     continue; 
                else
                     line(60:67) = spacereplace(line(60:67)); 
                     epoch_data(j,6) = str2num(line(60:67));
                end  
                line = fgetl(file);
            else
                line = fgetl(file);
            end
       
    end
    ind = find(sum(epoch_data,2)==0) ;
    epoch_data(ind,:) = [];

    curr_line = line;
    data{epoch,2} = epoch_data;
    epoch = epoch +1;
    n = n + 1;
end
obs.data = data;
disp('----------------Completed reading obs file-------------------');
fclose(file);
end