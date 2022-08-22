function [rcvalues] = countblackspots(matrix)
%In matrix 'eye', everything with a value of 1 is a black color.
%An eye should exist of more pixels then 1, so more pixels next to
%eachother with a value of 1 could be the eye.

%If the sum of 1 pixel plus it's surrounding pixels = 9, then this could be
%the eye.

[a b] = size(matrix);
rvalues = [];
cvalues = [];
for m = 2:a-1
    for n = 2:b-1
        if matrix(m,n)== 1 
          rvalues = [rvalues m];
          cvalues = [cvalues n];
        end   
    end
end
rcvalues = [rvalues;cvalues];
end

         