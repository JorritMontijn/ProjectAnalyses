function [x_vis,y_vis,vispar,x_au,y_au,aupar] = MOL_Psy_Trials(varargin)
%% Get input arguments:
tempsessionData     = varargin{1};
temptrialData       = varargin{2};

%% Get the trialtypes per session:
trialtypes = unique(temptrialData.trialType);

%% Define probe, auditory and visual trials and associated settings:
seshasprobe = 0; seshasau = 0; seshasvis = 0;
for tt = 1:length(trialtypes)
    switch trialtypes{tt}
        case 'P'
            seshasprobe         = 1;
            probetrialtype      = 'P';
        case 'Q'
            seshasprobe         = 1;
            probetrialtype      = 'Q';
        case 'R'
            seshasprobe         = 1;
            probetrialtype      = 'R';
        case 'V'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasvis           = 1;
                visuallog          	= 1;
                vistrialtype        = 'V';
                vismultiply100      = 1;
                visxaxislabel       = 'Contrast (%)';
                vissessionfield     = 'vecVisStimIntensities';
                vistrialfield       = 'visualInt';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasvis           = 1;
                vismultiply100      = 0;
                visuallog           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecOriChange';
                vistrialfield       = 'visualOriChange';
                visxaxislabel       = 'Change in orientation (degrees)';
            end
        case 'X'
            seshasvis           = 1;
            vismultiply100      = 0;
            visuallog           = 1;
            vistrialtype        = 'X';
            vissessionfield     = 'vecOriChange';
            vistrialfield       = 'visualOriChange';
            visxaxislabel       = 'Change in orientation (degrees)';
        case 'A'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasau            = 1;
                audiolog            = 0;
                aumultiply100       = 1;
                autrialtype         = 'A';
                ausessionfield      = 'vecAuStimIntensities';
                autrialfield        = 'audioInt';
                auxaxislabel        = 'Sound level (dB)';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasau            = 1;
                aumultiply100       = 0;
                audiolog            = 1;
                autrialtype         = 'A';
                ausessionfield      = 'vecFreqChange';
                autrialfield        = 'audioFreqChange';
                auxaxislabel        = 'Change in frequency (Hz)';
            end
        case 'Y'
            seshasau            = 1;
            aumultiply100       = 0;
            audiolog            = 1;
            autrialtype         = 'Y';
            ausessionfield      = 'vecFreqChange';
            autrialfield        = 'audioFreqChange';
            auxaxislabel        = 'Change in frequency (Hz)';
    end
end

if seshasvis
    if iscell(tempsessionData.(vissessionfield))
        visconditions   = tempsessionData.(vissessionfield){1};
    else
        visconditions   = tempsessionData.(vissessionfield);
    end
end

if seshasau
    if iscell(tempsessionData.(ausessionfield))
        auconditions    = tempsessionData.(ausessionfield){1};
    else
        auconditions    = tempsessionData.(ausessionfield);
    end
end

%% If session has opto:
selectionopto = [];
% if isfield(tempsessionData,'UseOpto')
%     seshasopto = tempsessionData.UseOpto;
% else         seshasopto = 0;
% end
seshasopto = 0;

if seshasopto && isfield(temptrialData,'hasphotostim') && any(temptrialData.hasphotostim==1)
    selectionopto(:,1) = ~(temptrialData.hasphotostim==1);
    selectionopto(:,2) = (temptrialData.hasphotostim==1);
else selectionopto = ones(length(temptrialData.trialStart),1);
end

%% Selection of responses per trial:
if tempsessionData.VisualLeftCorrectSide
    responseasauditory  = strcmp(temptrialData.responseSide,'R');
    responseasvisual    = strcmp(temptrialData.responseSide,'L');
else
    responseasauditory  = strcmp(temptrialData.responseSide,'L');
    responseasvisual    = strcmp(temptrialData.responseSide,'R');
end

%% Selection of probe trials:
if seshasprobe
    selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
    fprintf('%d Probe Trials\n',sum(selectionprobe))
else selectionprobe = strcmp(temptrialData.trialType,'Z'); %Nonsensical index
end

%% Auditory:
if seshasau
    for iOpto = 1:size(selectionopto,2)
        % Auditory trials:
        selectionau         = NaN(length(temptrialData.trialType),length(auconditions)+1);
        selectionau(:,1)    = selectionprobe & selectionopto(:,iOpto);
        for i = 1:length(auconditions)
            %                 selectionau(:,i+1)    = strcmp(temptrialData.trialType,autrialtype) & ismember(-temptrialData.(autrialfield),auconditions(i));
            selectionau(:,i+1)    = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(i)) & selectionopto(:,iOpto);
            fprintf('%d Audio Trials (%d) %s\n',sum(selectionau(:,i+1)),auconditions(i),repmat('Opto',1,iOpto-1))
        end
        
        percsidecorrectauditory = NaN(1,size(selectionau,2));
        percsideincorrectauditory = NaN(1,size(selectionau,2));
        for i = 1:size(selectionau,2)
            percsidecorrectauditory(i)    = sum(selectionau(:,i) & responseasauditory) / sum(selectionau(:,i));
            percsideincorrectauditory(i)    = sum(selectionau(:,i) & responseasvisual) / sum(selectionau(:,i));
        end
        
        if aumultiply100; auconditions = auconditions*100; end
        
        if numel(auconditions)==1; audiolog = 0; end
        
        if audiolog %If logarithmic axis place probe -2 times interval on log scale
            auprobepos = exp(min(log(auconditions)) - 2*mean(diff(log(auconditions))));
        else auprobepos = 0;
        end
        
        % Construct x from position of the probe trial and conditions
        xdata_au = [auprobepos auconditions];
        
        %% Psychometric fit:
        if numel(xdata_au)>=4 %If sufficient data points: fit a psychometric curve
            if audiolog
                xaulogdata = log10(xdata_au);
                [mu, sd, gr, lr, aucurve] = MOL_FitCumGauss(xaulogdata,percsidecorrectauditory);
                aucurve(:,1)  =   10.^aucurve(:,1);
                mu = 10.^mu;
                sd = 10.^sd;
            else
                [mu, sd, gr, lr, aucurve] = MOL_FitCumGauss(xdata_au,percsidecorrectauditory);
            end
            aupar.mu = mu;
            aupar.sd = sd;
            aupar.gr = gr;
            aupar.lr = lr;
            aupar.curve = aucurve;
        end
    end
end

%% Visual:
if seshasvis
    for iOpto = 1:size(selectionopto,2)
        selectionvis        = NaN(length(temptrialData.trialType),length(visconditions)+1);
        selectionvis(:,1)   = selectionprobe & selectionopto(:,iOpto);
        for i = 1:length(visconditions)
            %                 selectionvis(:,i+1)    = strcmp(temptrialData.trialType,vistrialtype) & ismember((temptrialData.(vistrialfield)),visconditions(i));
            selectionvis(:,i+1)    = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(i)) & selectionopto(:,iOpto);
            fprintf('%d Visual Trials (%d) %s\n',sum(selectionvis(:,i+1)),visconditions(i),repmat('Opto',1,iOpto-1))
        end
        percsidecorrectvisual   = NaN(1,size(selectionvis,2));
        percsideincorrectvisual = NaN(1,size(selectionvis,2));
        for i = 1:size(selectionvis,2)
            percsidecorrectvisual(i)    = sum(selectionvis(:,i) & responseasvisual ) / sum(selectionvis(:,i));
            percsideincorrectvisual(i)    = sum(selectionvis(:,i) & responseasauditory) / sum(selectionvis(:,i));
        end
        
        if vismultiply100; visconditions = visconditions*100; end
        
        if numel(visconditions)==1; visuallog = 0; end
        
        if visuallog %If logarithmic axis place probe -2 times interval on log scale
            visprobepos = exp(min(log(visconditions)) - 2*mean(diff(log(visconditions))));
        else visprobepos = 0;
        end
        
        % Construct x from position of the probe trial and conditions
        xdata_vis = [visprobepos visconditions];
        
        %% Psychometric fit:
        if numel(xdata_vis)>=4 %If sufficient data points: fit a psychometric curve
            if visuallog
                xvislogdata = log10(xdata_vis);
                [mu, sd, gr, lr, viscurve] = MOL_FitCumGauss(xvislogdata,percsidecorrectvisual);
                viscurve(:,1)  =   10.^viscurve(:,1);
                mu = 10.^mu;
                sd = 10.^sd;
            else
                [mu, sd, gr, lr, viscurve] = MOL_FitCumGauss(xdata_vis,percsidecorrectvisual);
            end
            vispar.mu = mu;
            vispar.sd = sd;
            vispar.gr = gr;
            vispar.lr = lr;
            vispar.curve = viscurve;
        end
    end
end


x_vis       = [0 visconditions];
y_vis       = percsidecorrectvisual;

x_au        = [0 auconditions];
y_au        = percsidecorrectauditory;

end