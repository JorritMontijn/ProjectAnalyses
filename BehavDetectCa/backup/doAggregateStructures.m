function sAggregate = doAggregateStructures(sAggregate,sTemp)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	cellFields = fieldnames(sTemp);
	
	for intField=1:length(cellFields)
		strField = cellFields{intField};
		if isfield(sAggregate,strField)
			if isstruct(sTemp.(strField))
				%sAggregate.(strField) = doAggregateStructures(sAggregate.(strField),sTemp.(strField));
			else
				sAggregate.(strField) = [sAggregate.(strField) sTemp.(strField)];
			end
		else
			sAggregate.(strField) = sTemp.(strField);
		end
	end
end

