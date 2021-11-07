matrix = [3 4 5 4; 6 0 4 5; 3 2 4 5; 33 0 4 5; 3 0 0 4];

size = size(matrix,1); 

for i = 1:size
    if(matrix(i,1) == 0)
        if(i == size)
            matrix([i-1,i],:) = matrix([i,i-1],:); 
            matrix(i-1,:) = []; 
            break;
        else
            matrix(i-1,:) = []; 
            i = i -1; 
            size = size - 1; 
        end
    elseif(matrix(i,2) == 0)
        if(i == size)
            matrix([i-1,i],:) = matrix([i,i-1],:); 
            matrix(i-1,:) = []; 
            break;
        else
            matrix(i,:) = []; 
            i = i -1; 
            size = size - 1; 
        end
%         matrix(i,:) = [];  
%         i = i -1; 
    elseif(matrix(i,3) == 0)
        if(i == size)
            matrix([i-1,i],:) = matrix([i,i-1],:); 
            matrix(i,:) = []; 
            break;
        else
            matrix(i,:) = []; 
            i = i -1; 
            size = size - 1; 
        end
%         matrix(i,:) = []; 
%         i = i -1; 
    elseif(matrix(i,4) == 0)
        if(i == size)
            matrix([i-1,i],:) = matrix([i,i-1],:); 
            matrix(i,:) = []; 
            break;
        else
            matrix(i-1,:) = []; 
            i = i -1; 
            size = size - 1; 
        end
%         matrix(i,:) = []; 
%         i = i -1;     
    else
        continue
    end
end