
function [plaidsact] = PlaidsactTuningData(orientationdata,contrastdata)
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



plaidsact.neuron   = sSesAggregate.neuron;
plaidsact.structStim   = sSesAggregate.structStim;

sSesAggregate.structStim.Orientation = sSesAggregate.structStim.Orientation(1:240);
sSesAggregate.structStim.TrialNumber = sSesAggregate.structStim.TrialNumber(1:240);
sSesAggregate.structStim.Contrast = sSesAggregate.structStim.Contrast(1:240);
sSesAggregate.structStim.FrameOn = sSesAggregate.structStim.FrameOn(1:240);
sSesAggregate.structStim.FrameOff = sSesAggregate.structStim.FrameOff(1:240);
%% Plaid laden

    plaid.trial = unique(sSesAggregate.structStim.TrialNumber(1:240));
    for i = 1:length(plaid.trial)
        plaid.contrindex(i).indices = find(sSesAggregate.structStim.TrialNumber(1:240) == plaid.trial(i));
        plaid.contrindex(i).orientation1 = sSesAggregate.structStim.Orientation(plaid.contrindex(i).indices(1));
        if length(plaid.contrindex(i).indices) == 2
            plaid.contrindex(i).orientation2 = sSesAggregate.structStim.Orientation(plaid.contrindex(i).indices(2));
        end
            plaid.contrindex(i).contrast1 = sSesAggregate.structStim.Contrast(plaid.contrindex(i).indices(1));
        if length(plaid.contrindex(i).indices) == 2
            plaid.contrindex(i).contrast2 = sSesAggregate.structStim.Contrast(plaid.contrindex(i).indices(2));
        end
    end
    for i = 1:length(plaid.contrindex)
      if (plaid.contrindex(i).orientation1 == 90)
        plaid.contrindex(i).orientation2 = 0;
      end
      if (plaid.contrindex(i).orientation1 == 0)
        plaid.contrindex(i).orientation2 = 90;
      end
      if isempty(plaid.contrindex(i).contrast2)
        plaid.contrindex(i).contrast2 = 0;
      end
    end
   for i = 1:length(plaid.contrindex)
        if length(plaid.contrindex(i).indices) == 2
            if sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(1)) > sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(2))
                plaid.contrindex(i).FrameOn = sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(1));
                plaid.contrindex(i).TrialfixOn = plaid.contrindex(i).indices(1);
            elseif sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(2)) >= sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(1))
                plaid.contrindex(i).FrameOn = sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices(2));
                plaid.contrindex(i).TrialfixOn = plaid.contrindex(i).indices(2);
            end
        elseif length(plaid.contrindex(i).indices) == 1
            plaid.contrindex(i).FrameOn = sSesAggregate.structStim.FrameOn(plaid.contrindex(i).indices);
            plaid.contrindex(i).TrialfixOn = plaid.contrindex(i).indices;
        end
            if length(plaid.contrindex(i).indices) == 2
            if sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(1)) < sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(2))
                plaid.contrindex(i).FrameOff = sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(1));
                plaid.contrindex(i).TrialfixOff = plaid.contrindex(i).indices(1);
            elseif sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(2)) <= sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(1))
                plaid.contrindex(i).FrameOff = sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices(2));
                plaid.contrindex(i).TrialfixOff = plaid.contrindex(i).indices(2);
            end
        elseif length(plaid.contrindex(i).indices) == 1
            plaid.contrindex(i).FrameOff = sSesAggregate.structStim.FrameOff(plaid.contrindex(i).indices);
            plaid.contrindex(i).TrialfixOff = plaid.contrindex(i).indices;
            end
   end
%%
plaidsact.UC = unique(plaidsact.structStim.Orientation); %vindt alle bestaande unieke oriëntaties uit de data
plaidsact.n_neuronen = length(vertcat(plaidsact.neuron.id));
plaidsact.n_contrasts = length(plaidsact.UC);
plaidsact.trials = plaid.trial;
plaidsact.n_trials = length(plaid.trial);


%preallocation for variables

plaidsact.activ = zeros(plaidsact.n_neuronen, plaidsact.n_trials);



    for i = 1:plaidsact.n_neuronen
        for k = 1:plaidsact.n_trials
            plaidsact.activ(i,k) = mean(plaidsact.neuron(i).dFoF(plaid.contrindex(k).FrameOn:plaid.contrindex(k).FrameOff));
        end
    end
plaidsact.activZ = zscore(plaidsact.activ,[],1); %z-scored at each trial to remove mean population activity
plaidsact.neuron   = sSesAggregate.neuron;
plaidsact.structStim   = sSesAggregate.structStim;
end

