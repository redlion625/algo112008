function obs = obs_read(filepath)
%D:\Third Year\ESSE 3670\Project 3\Data Downloaded\datasets to try with\algo\algo112008\ALGO0010.08O
file = fopen("Pixel4_derived.csv"); %"ALGO0010.08O"
disp('------------------Begin reading obs file---------------------');
%read constant parameters
curr_line = fgetl(file); %fget1 %fopen %fgetl
obs.rec_ID = []; 
obs.pos_xyz = [];
obs.ant_delta_HEN = [];
obs.num_typeobs = [];
obs.type_obs = [];
% 
% while isempty(strfind(curr_line,'END OF HEADER'))
%     if ~isempty(strfind(curr_line,'REC # / TYPE / VERS'))
%         temp = strsplit(curr_line);
%         obs.rec_ID = temp(1);
%     elseif ~isempty(strfind(curr_line,'APPROX POSITION XYZ'))
%         temp = strsplit(curr_line);
%         for j=1:3
%             obs.pos_xyz(j) = str2double([temp{j+1}]);
%         end
%     elseif ~isempty(strfind(curr_line,'ANTENNA: DELTA H/E/N'))
%         temp = strsplit(curr_line);
%         for j=1:3
%             obs.ant_delta_HEN(j) = str2double([temp{j+1}]);
%         end
%     elseif ~isempty(strfind(curr_line,'# / TYPES OF OBSERV'))
%         temp = strsplit(curr_line);
%         obs.num_typeobs = str2double(temp(2));
%         obs.type_obs = temp(3:obs.num_typeobs+2);
%     end
%     curr_line = fgetl(file); %fgetl or fget1
% end


disp('----------------Completed reading obs file-------------------');
fclose(file);
end