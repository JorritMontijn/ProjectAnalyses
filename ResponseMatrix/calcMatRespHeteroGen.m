function [vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matResp,sParams)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%get parameters
	if nargout == 2,boolDoAct = true;vecActivity=[];else boolDoAct = false;vecActivity=[];end
	if nargin < 2,sParams = struct;end
	
	%retrieve data from input
	intNeurons = size(matResp,1);
	intMaxT = size(matResp,2);
	
	%check if we should use vectorized or loop-based calculations 
	intDoubleBytes = 8;
	intMemReq = intDoubleBytes * intNeurons^2 * intMaxT * 3 * 2;
	[userview systemview] = memory; %#ok<ASGLU>
	if intMemReq > systemview.PhysicalMemory.Total
		warning([mfilename ':InsufMemForVecCalc'],'Insufficient memory for vectorized calculation; switching to loop-based calculation');
		intType = 1;
	else
		intType = 2;
	end
	
	%calculate heterogeneity
	if intType == 1 %slower, but requires less memory
		%perform non-vectorized calculation of heterogeneity over all time points
		matZ = squeeze(zscore(matResp,[],2)); %calculate vector-by-time
		if boolDoAct,vecActivity = squeeze(mean(matResp,1));end
		clear matResp; %clear unused variable
		vecHeterogeneity = nan(intMaxT,1);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		parfor intFrame=1:intMaxT
			%if mod(intFrame,100) == 0,fprintf('%s: Now at frame %d/%d\n',mfilename,intFrame,intMaxT);end
			vecZ = matZ(:,intFrame);
			matDist = abs(bsxfun(@minus,vecZ,vecZ'));
			vecHeterogeneity(intFrame) = mean(matDist(matSelect));
			if mod(intFrame,1000)==0,fprintf('now at %d/%d [%s]\n',intFrame,intMaxT,getTime);drawnow;end
		end
	else %faster, but requires more memory
		%perform vectorized calculation of heterogeneity over all time points
		matZ1 = nan(intNeurons,1,intMaxT); %create vertical vector-by-time
		matZ1(:,1,:) = zscore(matResp,[],2); %calculate vertical vector-by-time
		if boolDoAct,vecActivity = squeeze(mean(matResp,1));end
		clear matResp; %clear unused variable
		matZ2 = nan(1,intNeurons,intMaxT); %create horizontal vector-by-time
		matZ2(1,:,:) = squeeze(matZ1); %assign values to properly shaped matrix
		%matZ1 = repmat(matZ1,[1 intNeurons 1]); %repeat vertical matrix horizontally
		%matZ2 = repmat(matZ2,[intNeurons 1 1]); %repeat horizontal matrix vertically
		%matDist = abs(matZ1 - matZ2); %calculate z-score distance
		matDist = abs(bsxfun(@minus,matZ1,matZ2));
		clear('matZ1','matZ2'); %clear unused variables
		matDist(repmat(tril(true(intNeurons,intNeurons),0),[1 1 intMaxT])) = nan; %set auto-distance + lower triangle to nan
		vecHeterogeneity = squeeze(nanmean(reshape(matDist,[intNeurons^2 1 intMaxT]),1));%get output after reshaping
		clear matDist; %clear unused variable
	end
	vecActivity=vecActivity';
end
