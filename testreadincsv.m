% N = readtable('ground_truth.csv'); 
% N = readmatrix('ground_truth.csv'); 
% 
% M = table2array(N);

function output = referencecoordiantes(groundtruth)

M = groundtruth; %readmatrix('ground_truth.csv'); 

Time = M(:,3);
Lat = M(:,4); 
Long = M(:,5); 
Height = M(:,6); 

for i = 1:length(Time)
    TimeCorr(i,1) = Time(i)/1000; 
end

Combined = zeros(length(Time),4);
Combined(:,1) = Lat; 
Combined(:,2) = Long; 
Combined(:,3) = Height; 
Combined(:,4) = TimeCorr;

output = Combined;


