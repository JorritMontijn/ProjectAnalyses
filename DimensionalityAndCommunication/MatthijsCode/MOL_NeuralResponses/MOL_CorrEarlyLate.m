function MOL_CorrEarlyLate(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.zscore               = 1;

%% 

params.AlignOn          = 'stimChange';      %On which timestamp to align as t=0

% [splits,colors]         = MOL_GUI_GetSplits(trialData);             %Get splits on the basis of the GUI

params.showIndFig       = 1;

%% Main loop to get tensor:
neuroncounter               = 1;

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    events_ts = temptrialData.(params.AlignOn);

    %% construct normalized psth matrix:
    for iNeuron = 1:length(tempspikeData.ts)
        %% Get the spikes for this neuron:
        spikes_ts = tempspikeData.ts{iNeuron};
        
        %% Get histogram per cell
        [edges,hist_mat]  = calc_psth(events_ts,spikes_ts,params);

        params.twin_resp_start  = 0e6;
        params.twin_resp_stop   = 0.2e6;
        [resp_1,~] = calc_resp_from_psth(edges,hist_mat,params);
        
        params.twin_resp_start  = 0.2e6;
        params.twin_resp_stop   = 0.4e6;
        [resp_2,~] = calc_resp_from_psth(edges,hist_mat,params);
        
        if mean(resp_1)>1 || mean(resp_2)>1
            figure; plot(mean(hist_mat,1)); hold on; plot(mean(hist_mat(:,edges>0 & edges<0.6e6,1)));
            [dip,xl,xu, ifault, gcm, lcm, mn, mj] = HartigansDipTest(mean(hist_mat(:,edges>0 & edges<0.5e6,1)));
            text(1.5e3,1,num2str(dip),'FontSize',20);
        end
        
        if isempty(resp_1)
            keyboard
        end
        [r,p] = corr([resp_1 resp_2]);
        r_all(neuroncounter,1) = r(1,2);
        p_all(neuroncounter,1) = p(1,2);
        
                edges = -5:0.25:10;
        
        [S,Xedges,Yedges]           = histcounts2(resp_1,resp_2,edges,edges);
        S                           = imgaussfilt(S,0.6);
        [A(neuroncounter,:)]        = Fit2dGaussian(S, 0);

        r1_sm_all(neuroncounter,1) = mean(resp_1(temptrialData.visualOriChangeNorm==1));
        r2_sm_all(neuroncounter,1) = mean(resp_2(temptrialData.visualOriChangeNorm==1));
        r1_big_all(neuroncounter,1) = mean(resp_1(temptrialData.visualOriChangeNorm==2));
        r2_big_all(neuroncounter,1) = mean(resp_2(temptrialData.visualOriChangeNorm==2));
        
        if p(1,2)<0.05
            if params.showIndFig > 2
                figure;
                idx = temptrialData.visualOriPostNorm<=2;
                scatter(resp_1(iNeuron,idx),resp_2(iNeuron,idx),'g','filled'); hold on;
                idx = temptrialData.visualOriPostNorm>2;
                scatter(resp_1(iNeuron,idx),resp_2(iNeuron,idx),'r','filled'); hold on;
                text(0,0,num2str(r(1,2)))
            end
        end
        neuroncounter = neuroncounter +1;
    end
end


%% 

% 2. 2D Rotated Gaussian function ( A requires 6 coefs ).
f = @(A,X) A(1)*exp( -(...
    ( X(:,:,1)*cos(A(6))-X(:,:,2)*sin(A(6)) - A(2)*cos(A(6))+A(4)*sin(A(6)) ).^2/(2*A(3)^2) + ... 
    ( X(:,:,1)*sin(A(6))+X(:,:,2)*cos(A(6)) - A(2)*sin(A(6))-A(4)*cos(A(6)) ).^2/(2*A(5)^2) ) );

[m,n] = size(S);
[x,y]=meshgrid(1:n,1:m); X=zeros(m,n,2); X(:,:,1)=x; X(:,:,2)=y;

figure; 
hold all;

for iNeuron = 1:length(r_all)
    if p_all(iNeuron,1)<0.01
        contour(Xedges(1:end-1),Yedges(1:end-1),f(A(iNeuron,:),X),1,'k-');
    end
%     patch
end


%% 

figure;
scatter(r1_sm_all,r2_sm_all,'b','filled'); hold on;
scatter(r1_big_all,r2_big_all,'r','filled'); hold on;



% plot lines
% hold on; plot(A(2),A(4),'+b',vx_h,vy_h,'.r',vx_v,vy_v,'.g');
% plot contour:


%     snakemat_splits{iSplit} = snakemat;
% end

end