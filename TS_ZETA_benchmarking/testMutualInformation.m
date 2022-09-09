%test mutual information
vecX =normrnd(0,1,1,10);
vecY = vecX;%0.5*normrnd(0,1,1,numel(vecX)); %normalized mutual information should be 0.5

%perform kde
n=min(2^ceil(log2((numel(vecX)+numel(vecY))/2+1)-1),2^10);
[bandwidth,density,X,Y]=kde2d(cat(1,vecX,vecY)',n);
density = density ./ sum(density(:));

matJointP = density ./ sum(density(:));
vecMarginalPY = sum(matJointP,1);
vecMarginalPX = sum(matJointP,2);
indRemY = vecMarginalPY<(1/max(n,numel(vecY)));
indRemX = vecMarginalPX<(1/max(n,numel(vecX)));
vecMarginalPY(indRemY) = [];
vecMarginalPX(indRemX) = [];
matJointP(indRemX,:) = [];
matJointP(:,indRemY) = [];

%renormalize
vecMarginalPY = vecMarginalPY ./sum(vecMarginalPY(:));
vecMarginalPX = vecMarginalPX ./sum(vecMarginalPX(:));
matJointP = matJointP ./sum(matJointP(:));

EntropyX = -sum(vecMarginalPX.*log2(vecMarginalPX));
EntropyY = -sum(vecMarginalPY.*log2(vecMarginalPY));

JointEntropyXY = -sum(sum(matJointP.*log2(matJointP)));

%calc MI
MutualInformation = EntropyX + EntropyY - JointEntropyXY;
NormalizedMutualInformation = (2*MutualInformation) / (EntropyX + EntropyY); %Witten & Frank, 2005

%run shuffles


%fit with beta distribution
boolFitWithBeta = true;
if ~(v<(mu*(1-mu)))
	disp('closed form parametrization with mean and variance is invalid');
	boolFitWithBeta = false;
end
if boolFitWithBeta
	%get parameters
	alpha = mu*((mu*(1-mu))/v - 1);
	beta = (1-mu)*((mu*(1-mu))/v - 1);
	%calculate p
	p = 1-betacdf(NormalizedMutualInformation,alpha,beta,'upper');
else
	%use quantiles
	
end

