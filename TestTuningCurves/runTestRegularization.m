


%% test script for regression
intNumQ = 10;
intT = 300;
vecQ1 = zeros(1,intT);
vecQ2 = zeros(1,intT);
vecQ2(101:200) = 1;
vecQ3 = zeros(1,intT);
vecQ3(201:intT) = 2;
vecRealB = vecQ1 + vecQ2 + vecQ3;
matY = nan(intNumQ,intT);
for intN=1:intNumQ
	matY(intN,:) = (vecRealB.*rand(size(vecRealB))+rand(size(vecRealB)));
end
matX = cat(1,vecQ1,vecQ2,vecQ3)';
matY = matY';

%%
dblLambda = 1;
[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matX,matY,1,dblLambda)
 
[matY_hat2,dblR2_CV2,matB2] = doCrossValidatedDecodingRR2(matX,matY,1,dblLambda)


%%

