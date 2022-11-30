%make some distributions
h=figure;maxfig;
vNum = 1000;
for intDistro=1:6
	%generate data
	if intDistro == 1
		%linear inverted
		
		vecX =normrnd(0,1,1,vNum);
		vecY =-vecX;%normrnd(0,1,1,vNum);%0.5*vecX + 0.5*normrnd(0,1,1,vNum); %normalized mutual information should be 0.5
		
	elseif intDistro == 2
		%X 1/10th scale of Y
		vecX =normrnd(0,1,1,vNum);
		vecY =0.1*vecX;%normrnd(0,1,1,vNum);%0.5*vecX + 0.5*normrnd(0,1,1,vNum); %normalized mutual information should be 0.5
		
	elseif intDistro == 3
		%uncorrelated
		vecX =normrnd(0,1,1,vNum);
		vecY =normrnd(0,1,1,vNum);%0.5*vecX + 0.5*normrnd(0,1,1,vNum); %normalized mutual information should be 0.5
		
	elseif intDistro == 4
		%X cross
		vecX =normrnd(0,1,1,vNum);
		vecSign = 2*double(rand(size(vecX))<0.5)-1;
		vecY =vecX.*vecSign;%normrnd(0,1,1,vNum);%0.5*vecX + 0.5*normrnd(0,1,1,vNum); %normalized mutual information should be 0.5
		
	elseif intDistro == 5
		%50 correlation
		vecX =normrnd(0,1,1,vNum);
		vecY =0.5*vecX + 0.5*normrnd(0,1,1,vNum);
	elseif intDistro == 6
		%sine wave
		vecX =normrnd(0,1,1,vNum);
		vecY =sin(vecX*10);
	end
	
	%perform kde
	vecRange = [-4 4];
	n=min(2^ceil(log2((vNum+vNum)/2+1)-1),2^10);
	[bandwidth,density,X,Y]=kde2d(cat(1,vecX,vecY)',n,vecRange(1).*[1 1],vecRange(2).*[1 1]);
	density = density + min(density(density>0));
	matJointP = density ./ sum(density(:));
	
	if 0
		%plot data
		figure;maxfig;
		subplot(2,3,1);
		dblStep = 0.25;
		vecRange = [-4 4];
		vecBinEdges = vecRange(1):dblStep:vecRange(2);
		vecBinCenters = vecBinEdges(2:end)-dblStep/2;
		cX = histcounts(vecX,vecBinEdges);
		plot(vecBinCenters,cX);
		title(sprintf('X, n=%d, normally distributed',vNum));
		xlabel('X value');
		ylabel('Count');
		fixfig;
		
		subplot(2,3,2);
		cY = histcounts(vecY,vecBinEdges);
		plot(vecBinCenters,cY);
		xlabel('Y value');
		ylabel('Count');title(sprintf('Y, n=%d, normally distributed',vNum));
		fixfig;
		
		subplot(2,3,3);
		scatter(vecX,vecY,100,'.');
		xlim(vecRange);
		ylim(vecRange);
		title(sprintf('Joint X,Y, Pearson r=%.3f',corr(vecX',vecY')));
		xlabel('X value');
		ylabel('Y value');
		fixfig;
		
		subplot(2,3,4);
		scatter(vecX(randperm(vNum)),vecY(randperm(vNum)),100,'.');
		xlim(vecRange);
		ylim(vecRange);
		title(sprintf('P(X)*P(Y)'));
		xlabel('X value');
		ylabel('Y value');
		fixfig;
		
		
		%perform stupid kde
		[f,xi,bw] = ksdensity(cat(1,vecX,vecY)');
		%kdeBandwidthX = 0.9*min(std(vecX),iqr(vecX)/1.34)*(vNum^-(1/5));
		%kdeBandwidthY = 0.9*min(std(vecX),iqr(vecY)/1.34)*(vNum^-(1/5));
		kdeBandwidthX = bw(1);
		kdeBandwidthY = bw(2);
		dblStepX = mean(median(diff(X,1,2)));
		dblStepY = mean(median(diff(Y,1,1)));
		gridKdeX = (-kdeBandwidthX*2):dblStepX:(kdeBandwidthX*2);
		kernelX = normpdf(gridKdeX,0,kdeBandwidthX);
		gridKdeY = (-kdeBandwidthY*2):dblStepY:(kdeBandwidthY*2);
		kernelY = normpdf(gridKdeY,0,kdeBandwidthY);
		kernel = kernelY' * kernelX;
		kernel = kernel ./ sum(kernel(:));
		
		gridX = X(1,:);
		gridY = Y(:,1)';
		gridEdgesX = [gridX(1) - dblStepX gridX] +dblStepX/2;
		gridEdgesY = [gridY(1) - dblStepY gridY] +dblStepY/2;
		
		[histGrid2,Xedges,Yedges] = histcounts2(vecX,vecY,gridEdgesX,gridEdgesY);
		stupidKde = imfilt(histGrid2,kernel);
		
		subplot(2,3,5)
		imagesc(vecRange,vecRange,stupidKde);
		title(sprintf('KDE with matlab''s ksdensity()'));
		xlabel('X value');
		ylabel('Y value');
		fixfig;grid off;
		
		subplot(2,3,6)
		imagesc(vecRange,vecRange,density);
		title(sprintf('Bivariate KDE using diffusion estimation'));
		xlabel('X value');
		ylabel('Y value');
		fixfig;grid off;
	end
	
	%calculate marginals
	vecMarginalPY = sum(matJointP,1);
	vecMarginalPX = sum(matJointP,2);
	indRemY = vecMarginalPY<(1/max(n,vNum));
	indRemX = vecMarginalPX<(1/max(n,vNum));
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
	%if ~(variance<(mu*(1-mu)))
	%	disp('closed form parametrization with mean and variance is invalid');
	%	boolFitWithBeta = false;
	%end
	if boolFitWithBeta
		%get parameters
		%alpha = mu*((mu*(1-mu))/v - 1);
		%beta = (1-mu)*((mu*(1-mu))/v - 1);
		%%calculate p
		%p = 1-betacdf(NormalizedMutualInformation,alpha,beta,'upper');
	else
		%use quantiles
		
	end
	
	
	% P(X ? x,Y ? y)
	%% https://epubs.siam.org/doi/epdf/10.1137/1.9781611972832.22
	%% https://doi.org/10.1007/s10618-022-00847-y
	%vecX =normrnd(0,1,1,1000);
	%vecY =normrnd(0,1,1,numel(vecX));%0.5*vecX + 0.5*normrnd(0,1,1,xNum); %normalized mutual information should be 0.5
	%vecY=vecX;
	vNum = numel(vecX);
	[vecXs,vecReorderX] = sort(vecX);
	vecInvertX = accumarray(vecReorderX',1:vNum);
	[vecYs,vecReorderY] = sort(vecY);
	vecInvertY = accumarray(vecReorderY',1:vNum);
	
	%marginals
	cdfV = linspace(1/vNum,1,vNum);
	cdfX = cdfV(vecInvertX);
	cdfY = cdfV(vecInvertY);
	
	%r=1, sum(joint)=3.85
	%r=-1, sum(joint)=2.2
	
	
	%cumulative entropy
	dX =[diff(vecXs) 0];
	Hx1 = -sum(...
		dX.*cdfV.*log(cdfV)...
		);
	
	dY = [diff(vecYs) 0];
	Hy1 = -sum(...
		dY.*cdfV.*log(cdfV)...
		);
	
	%error('test joint pdf: is joint cdf cumsum(cumsum(jointPdf,1),2)?')
	%conditional entropy
	cdfJoint = nan(vNum,vNum);
	H_full_x = nan(vNum,vNum);
	H_full_y = nan(vNum,vNum);
	H_full_yx = nan(vNum,vNum);
	H_full_xy = nan(vNum,vNum);
	Py_given_x = nan(vNum,vNum);
	Px_given_y = nan(vNum,vNum);
	CMI_full = nan(vNum,vNum);
	CMI_numer = nan(vNum,vNum);
	CMI_denom = nan(vNum,vNum);
	
	%normalize
	for yi=1:vNum
		for xi=1:vNum
			cdfJoint(xi,yi) = sum(vecX<=vecXs(xi) & vecY<=vecYs(yi))./vNum;
			Py_given_x(xi,yi) = cdfJoint(xi,yi) / (sum(vecXs<=vecXs(xi))./vNum);
			Px_given_y(xi,yi) = cdfJoint(xi,yi) / (sum(vecYs<=vecYs(yi))./vNum);
			
			Py = sum(vecYs<=vecYs(yi))./vNum;
			H_full_x(xi,yi) = dY(yi)*dX(xi)*cdfJoint(xi,yi)*log(Py);
			H_full_y(xi,yi) = dY(yi)*dX(xi)*cdfJoint(xi,yi)*log(Py);
			H_full_yx(xi,yi) = dY(yi)*dX(xi)*cdfJoint(xi,yi)*log(Py_given_x(xi,yi));
			H_full_xy(xi,yi) = dY(yi)*dX(xi)*cdfJoint(xi,yi)*log(Px_given_y(xi,yi));
			
			CMI_numer(xi,yi) = dX(xi)*dY(yi)*cdfJoint(xi,yi)*log(Py_given_x(xi,yi));
			CMI_denom(xi,yi) = dX(xi)*dY(yi)*Py*log(Py);
			
			
		end
	end
	CMI_numer(isnan(CMI_numer))=0;
	CMI_denom(isnan(CMI_denom))=0;
	Hx = -sum(sum(H_full_x));
	Hy = -sum(sum(H_full_y));
	Hxy = -sum(sum(H_full_xy));
	Hyx = -sum(sum(H_full_yx));
	H_joint = Hyx + Hx;
	CMI = 1 - mean(sum(CMI_numer(1:(end-1)),2)./sum(CMI_denom(1:(end-1)),2));
	
	if 0
		figure
		subplot(2,3,1)
		scatter(vecX,cdfX);
		title('Marginal cdf P(x)');
		fixfig;
		xlabel('X');
		ylabel('cumulative P(X)');
		fixfig;
		
		subplot(2,3,2)
		scatter(vecY,cdfY);
		title('Marginal cdf P(y)');
		xlabel('Y');
		ylabel('cumulative P(Y)');
		fixfig;
		
		subplot(2,3,3)
		scatter(vecX,vecY);
		title('Sample x,y');
		xlabel('X');
		ylabel('Y');
		fixfig;
		
		subplot(2,3,4)
		[meshX,meshY]=meshgrid(vecXs,vecYs);
		surf(meshX,meshY,Px_given_y,'EdgeColor','none');
		axis xy;
		view(90,90)
		title('Conditional cdf P(X<=x|y)');
		xlabel('X');
		ylabel('X');
		fixfig;grid off;
		
		subplot(2,3,5)
		surf(meshX,meshY,Py_given_x,'EdgeColor','none');
		axis xy;
		view(90,90)
		title('Conditional cdf P(Y<=y|x)');
		xlabel('X');
		ylabel('Y');
		fixfig;grid off;
		
		pdfJoint = diff(cdfJoint(:,2:end),1,1) + diff(cdfJoint(2:end,:),1,2);
		pdfJoint = cat(2,zeros(size(pdfJoint,1)+1,1),cat(1,zeros(1,size(pdfJoint,1)),pdfJoint));
		subplot(2,3,6)
		surf(meshX,meshY,log(cdfJoint),'EdgeColor','none');
		axis xy;
		view(90,90)
		title(sprintf('Joint cdf P(x,y); CMI=%.3f',CMI));
		xlabel('X');
		ylabel('Y');
		fixfig;grid off;
	end
	
	%% make plot
	subplot(2,3,intDistro,'Parent',h);
	scatter(vecX,vecY,'.');
	fixfig;
	xlabel('X');
	ylabel('Y');
	xlim([-4 4]);
	ylim([-4 4]);
	title(sprintf('Pearson r=%.3f, NMI=%.3f, CMI=%.3f',corr(vecX',vecY'),NormalizedMutualInformation,CMI));
end