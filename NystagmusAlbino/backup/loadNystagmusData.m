%loadNystagmusData Loads data from nystagmus project using a set of
%requested parameters (see defineNystagmusData.m), including:
%
%strSelectMouseType		'WT' or 'Albino'
%strSelectArea			'V1', 'NOT', or 'SC'
%strSelectMouse			e.g., 'MB2' or 'MA1'
%strSelectDate			e.g., '20190315'
%strSelectStim			'RF', 'DG', or 'RF-DG'
%intSelectPopulation	Integer specifying unique neuronal population #
%
%Version History:
%2019-06-11 Created [by Jorrit Montijn]

%% define data
defineNystagmusData;

%% select sessions
indSelectSession = true(size(sSession));
cellSelectBy = {'strSelectMouseType','strSelectArea','strSelectMouse','strSelectDate'};
for intSelectBy=1:numel(cellSelectBy)
	strQueryVar = cellSelectBy{intSelectBy};
	strProperty = strQueryVar(10:end);
	strQueryVar = strcat('strSelect',strProperty);
	strEvalLeft = strcat('ind',strProperty,'=');
	if exist(strQueryVar,'var') && ~isempty(eval(strQueryVar))
		%build code
		strEvalRight = strcat('strcmpi({sSession(:).',strProperty,'},',strQueryVar,');');
		%assign to index vector
		eval(strcat(strEvalLeft,strEvalRight));
		%narrow selection index
		eval(strcat('indSelectSession = indSelectSession & ',strEvalRight));
	end
end
%% select recordings
%strSelectStim = strStim;
%intSelectPopulation = intPopulation;
cellFields = {'Location','Depth'};
intPopCounter = 0;
boolFound = false;
boolNonUnique = false;
for intSession=find(indSelectSession)
	%get recording structure
	sRec = sSession(intSession).Recording;
	%get unique types
	sUniqueCombos = getUniqueFieldCombos(sRec,{'Location','Depth'});
	vecNeuralPopulation = sUniqueCombos.vecUniqueIdx;
	
	%get requested stimulus
	indReqStim = ~cellfun(@isempty,strfind({sRec.Stim},strSelectStim)); %#ok<STRCLFH>
	if any(indReqStim)
		vecThesePops = vecNeuralPopulation(indReqStim);
		vecIndexedPops = label2idx(vecThesePops) + intPopCounter;
		vecReqStim = find(indReqStim);
		vecSelectRec = vecReqStim(vecIndexedPops(:) == intSelectPopulation);
		
		%increment counter
		intPopCounter = intPopCounter + numel(vecIndexedPops);
		
		%check if found
		if numel(vecSelectRec) > 0
			if numel(vecSelectRec) == 1 && ~boolFound
				%get data
				boolFound = true;
				strBlock = sRec(vecSelectRec).Block;
				strMouse = sSession(intSession).Mouse;
				strDate = sSession(intSession).Date;
				intSelectSession = intSession;
				intSelectRec = vecSelectRec;
			else
				boolNonUnique = true;
			end
		end
	end
end


%% error messages
if ~boolFound
	error([mfilename ':FileNotFound'],sprintf('Recording (property combination) not found, max pop=%d',intPopCounter));
elseif boolNonUnique
	error([mfilename ':MultipleResults'],sprintf('Non-unique combination of properties requested, max pop=%d',intPopCounter));
end

%% set paths
strBlock = num2str(strBlock);
strDataRoot = 'D:\Data\Raw\ePhys';
strStimLog = [strDataRoot filesep 'StimLogs' filesep strMouse '_' strDate];

% set paths
ops.root = [strDataRoot filesep 'KiloSortBinaries']; % 'openEphys' only: where raw files are
ops.rec  = [strMouse '_' strDate '_B' strBlock]; %which recording to process


