function [contrastact] = ContrastactTuningData(orientationdata,contrastdata)
%% concatenate
%add contrast field to gratings
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

contrastact.neuron   = sSesAggregate.neuron;
contrastact.structStim   = sSesAggregate.structStim;

contrastact.UC = unique(contrastact.structStim.Contrast); %vindt alle bestaande unieke oriëntaties uit de data
contrastact.n_neuronen = length(vertcat(contrastact.neuron.id));
contrastact.n_contrasts = length(contrastact.UC);
contrastact.trials = contrastact.structStim.Contrast;
contrastact.n_trials = length(contrastact.trials);


%preallocation for variables


contrastact.stimindex = struct("indices", cell(1, contrastact.n_contrasts), 'contrast', cell(1,contrastact.n_contrasts));
contrastact.neur = struct('meansvect', cell(1, contrastact.n_contrasts));
contrastact.neurvect = zeros(1,contrastact.n_neuronen);
contrastact.max_val = zeros(1,contrastact.n_trials);
contrastact.activ = zeros(contrastact.n_neuronen, contrastact.n_trials);

%finding the different stimulus contrasts and their respective indexes.

    for i = 1:contrastact.n_contrasts
        contrastact.stimindex(i).indices=find(contrastact.structStim.Contrast==contrastact.UC(i));
        contrastact.stimindex(i).contrast=contrastact.UC(i);
    end


%determining the mean response of every neuron to the different stimulus
%presentations, corresponding to the contrast of that stimulus
    for k = 1:contrastact.n_neuronen
        for i = 1:contrastact.n_contrasts
            for j = 1:length(contrastact.stimindex(i).indices)
                contrastact.stimindex(i).neuron(k).meanresponse(j) = mean(contrastact.neuron(k).dFoF(contrastact.structStim.FrameOn(contrastact.stimindex(i).indices(j)):contrastact.structStim.FrameOff(contrastact.stimindex(i).indices(j))));
            end
        end
    end
  % voor elke oriëntatie worden voor alle neuronen de gemiddelde activiteiten tijdens stimuluspresentatie bepaald en als .meanresponse opgeslagen.


    for k = 1:contrastact.n_neuronen
        for i = 1:contrastact.n_contrasts
            contrastact.stimindex(i).neuron(k).grandmean = mean(contrastact.stimindex(i).neuron(k).meanresponse);%voor alle contrasten wordt voor alle neuronen het gemiddelde berekend van de gemiddelde activiteiten tijdens stimuluspresentatie.
            contrastact.neur(k).meansvect(i) = contrastact.stimindex(i).neuron(k).grandmean; %maakt een nieuwe struct aan voor de gemiddeldes. Maakt vorige stap arbitrair?
            contrastact.neurvect(i,k) = contrastact.neur(k).meansvect(i); %maakt een matrix met de gemiddeldes, waarbij de rows de orientaties zijn en de collumns de neuronen
        end
    end



    for i = 1:contrastact.n_neuronen
        for k = 1:contrastact.n_trials
            contrastact.activ(i,k) = mean(contrastact.neuron(i).dFoF(contrastact.structStim.FrameOn(k):contrastact.structStim.FrameOff((k))));
        end
    end
contrastact.activZ = zscore(contrastact.activ,[],2); %z-scored at each trial to remove mean population activity


%matresp = getRespMat(sSesAggregate);





end

