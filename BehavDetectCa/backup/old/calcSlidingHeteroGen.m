function [vecHeterogeneity,vecActivity] = calcSlidingHeteroGen(ses,sParams)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%get parameters
	if nargout == 2,boolDoAct = true;vecActivity=[];else boolDoAct = false;vecActivity=[];end
	if nargin < 2,sParams = struct;end
	if isfield(sParams,'intWindowLength'),intWindowLength=round(sParams.intWindowLength);else intWindowLength = 1;end
	if isfield(sParams,'intWindowLength') && intWindowLength ~= sParams.intWindowLength, warning([mfilename ':InputNotAnInteger'],'sParams.intWindowLength was supplied as non-integer value [%.3f]; will use [%d] instead',sParams.intWindowLength,intWindowLength);end
	
	%retrieve data from input structure and assign to response matrix
	intNeurons = numel(ses.neuron);
	intMaxT = length(ses.neuron(1).dFoF);
	matResp = nan(intNeurons,1,intMaxT);
	vecKernel = ones(1,intWindowLength)/intWindowLength;
	for intNeuron=1:intNeurons
		if intWindowLength > 1 %perform sliding window averaging if required
			matResp(intNeuron,1,:) = [nan(1,floor(intWindowLength/2)) conv(ses.neuron(intNeuron).dFoF,vecKernel,'valid') nan(1,floor(intWindowLength/2))];
		else
			matResp(intNeuron,1,:) = ses.neuron(intNeuron).dFoF;
		end
	end
	
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
		matZ = squeeze(zscore(matResp,[],3)); %calculate vector-by-time
		if boolDoAct,vecActivity = squeeze(mean(matResp,1));end
		clear matResp; %clear unused variable
		vecHeterogeneity = nan(intMaxT,1);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		for intFrame=1:intMaxT
			vecZ = matZ(:,intFrame);
			matDist = abs(bsxfun(@minus,vecZ,vecZ'));
			vecHeterogeneity(intFrame) = mean(matDist(matSelect));
		end
	else %faster, but requires more memory
		%perform vectorized calculation of heterogeneity over all time points
		matZ1 = zscore(matResp,[],3); %calculate vertical vector-by-time
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
end

