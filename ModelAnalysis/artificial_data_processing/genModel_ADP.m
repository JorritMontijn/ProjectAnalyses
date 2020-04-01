%function varOut = genModel_ADP(strType,matParams)
strType = 'poly';
matParams = (1:10)';
x=1:10;

	cellModelTypes = {'poly','gauss','sigmoid'};
	
	if strcmpi(strType,cellModelTypes{1}) %poly
		varBasisFunction = @(x,w,matParams) bsxfun(@times,w,bsxfun(@power,x,matParams(:,1)-1));
	elseif strcmpi(strType,cellModelTypes{2}) %gauss
		varBasisFunction = @(x,w,matParams) exp( -(x - matParams(1,:)).^2 ./ (2*matParams(2,:).^2)); %matParams: {mu, s}
	elseif strcmpi(strType,cellModelTypes{3}) %sigmoid
		varBasisFunction = @(x,w,matParams) (1 / (1 + exp(-matParams(3,:))) ) .* ((x - matParams(1,:)) ./ matParams(2,:)); %matParams: {mu, s, a}
	else
		
	end

%end
