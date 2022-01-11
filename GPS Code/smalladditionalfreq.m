takein = dlmread('testmatrix.txt');
obssize = size(takein,1);
f1 = 1575.42e6; %L1
f2 = 1176.45e6; %L2
 for j = 1:obssize(1) %Iterating through the submatrix attached to the observation epoch
         %P2 = Observations{i,2}(j,5); %P1 code
         %P1 = Observations{i,2}(j,6); %P2 code
         P1 = takein(j,3); 
         P2 = takein(j,7);
         if(P2 == 0)
             Pseudorange(j) = P1; %Observations{i,2}(j,3);
         else
             Pseudorange(j) = ((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2);  %P1;%((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2);
             testPseudorange(j) =  ((f1^2)/((f1^2) - (f2^2)))*(P1) - ((f2^2)/((f1^2) - (f2^2)))*(P2);
         end
 end
 
 for i = 1:obssize(1)
     result(i) = Pseudorange(i) - testPseudorange(i); 
 end
 
 testpseudo = P1; 
 testdiff = Pseudorange(6) - P1;
%          Pseudorange(j) = Observations{i,2}(j,3); %((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2); %Calculating the uncorrected pseudorange
%          index = Observations{i,2}(j,11); 
%          sqrtA = Nav_data.data(index, 13); %Finding the semi-major axis
%          A = (sqrtA)^2; %Finding the semi-major axis
%          N0 = sqrt(gravitationalparameter/A^3); %Computing the mean motion
%  end