function [orientact] = OrientactTuningData(orientationdata,contrastdata)
%% concatenate
%add contrast field to gratings
sesG.structStim.Contrast = 100*ones(size(orientationdata.ses.structStim.Orientation));
sSesAggregate = buildMultiSesAggregate(orientationdata.ses, contrastdata.ses);

for intCellType=1:numel(sSesAggregate.cellType)
	strCellType = sSesAggregate.cellType{intCellType};
	if isfield(orientationdata.ses,strCellType)
		indKeepG = ~cellfun(@strcmpi,cellfill('absent',[1 numel(orientationdata.ses.(strCellType))]),{orientationdata.ses.(strCellType).strPresence});
		indKeepP = ~cellfun(@strcmpi,cellfill('absent',[1 numel(contrastdata.ses.(strCellType))]),{contrastdata.ses.(strCellType).strPresence});
		
		%check which to keep
		vecKeepAgg = zeros([1 numel(sSesAggregate.(strCellType))]);
		vecKeepAgg(indKeepG) = 1;
		vecKeepAgg(indKeepP) = vecKeepAgg(indKeepP) + 1;
		indKeepAgg = vecKeepAgg==2;
		%remove entries
		sSesAggregate.(strCellType) = sSesAggregate.(strCellType)(indKeepAgg);
	end
end

orientact.neuron   = sSesAggregate.neuron;
orientact.structStim   = sSesAggregate.structStim;
orientact.n_neuronen = length(vertcat(orientact.neuron.id));
orientact.trials = orientact.structStim.Orientation;
orientact.n_trials = length(orientact.trials);

%preallocation for variables

activa = zeros(orientact.n_neuronen, orientact.n_trials);


for i = 1:orientact.n_neuronen
    for k = 1:orientact.n_trials
        activa(i,k) = mean(orientact.neuron(i).dFoF(orientact.structStim.FrameOn((k)):orientact.structStim.FrameOff((k))));
    end
end


indGratings = [241:320];

orientact.activ = activa(:,indGratings);
orientact.activZG = zscore(orientact.activ, [], 2);

end


