function matAggResp = getPlaidData(sRec,strComparison,boolPV)
	%getPlaidData Summary of this function goes here
	%    matAggResp = getPlaidData(sRec,strComparison,boolPV)
	
	if ~exist('boolPV','var') || isempty(boolPV)
		boolPV = false;
	end
	if boolPV
		strExpHeader = 'matRespPV(:,';
	else
		strExpHeader = 'matResp(:,';
	end
	
	matAggResp = [];
	for intRec=1:numel(sRec)
		sThisRec = sRec(intRec);
		%unpack structure
		cellFields = fieldnames(sThisRec);
		for intField=1:numel(cellFields)
			strField = cellFields{intField};
			eval([strField '=sThisRec.' strField ';']);
		end
		
		if boolPV && isempty(matRespPV)
			continue;
		end
		%get data
		matAggResp = cat(1,matAggResp,eval([strExpHeader strComparison ');']));
	end
end

