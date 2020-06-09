function [OSI,DSI,gOSI] = calc_OSIDSI(all_oris,ori_resp,showFig)
%calculates orientation selectivity index

%Gets as input:
% 1 The orientations in degrees
% 2 The responses of an individual cell to each orientation
% 3 Boolean flag whether to plot figure
%Note that inputs 1 and 2 are of same size (vector 1xN)

%% Find preferred and opposite+orthogonal responses:
nOris           = length(all_oris);
ori_pref        = find(ori_resp==max(ori_resp)); ori_pref = ori_pref(1);
ori_ortho1      = mod(ori_pref+nOris/4,nOris); if ori_ortho1==0;  ori_ortho1 = nOris; end
ori_ortho2      = mod(ori_pref-nOris/4,nOris); if ori_ortho2==0;  ori_ortho2 = nOris; end
ori_oppo        = mod(ori_pref+nOris/2,nOris); if ori_oppo==0;    ori_oppo = nOris;   end

%% Caluclate OSI
OSI = (ori_resp(ori_pref)-mean([ori_resp(ori_ortho1) ori_resp(ori_ortho2)])) / ...
    (ori_resp(ori_pref)+mean([ori_resp(ori_ortho1) ori_resp(ori_ortho2)]));

%% Calculate DSI
DSI = (ori_resp(ori_pref)-ori_resp(ori_oppo)) / ...
    (ori_resp(ori_pref)+ori_resp(ori_oppo));

%% gOSI: 
gOSI = 1-circ_var((pi/180)*all_oris, ori_resp);

%this same term is applied to gDSI:
% gDSI: 
% gDSI = %normalized vector sum
% Shi et al 2017 - retinal origin of direction selectivity
%Which is mean resultant 

%% Plot if requested
if showFig
    figure; set(gcf,'color','w');
    plot(all_oris,ori_resp,'k.','MarkerSize',15); hold on;
    
    plot([all_oris(ori_pref) all_oris(ori_ortho1)],[max(ori_resp) mean([ori_resp(ori_ortho1) ori_resp(ori_ortho2)])],'b:','LineWidth',2); hold on;
    plot([all_oris(ori_pref) all_oris(ori_oppo)],[max(ori_resp) ori_resp(ori_oppo)],'r:','LineWidth',2); hold on;
    title(sprintf('OSI is %2.2f, DSI is %2.2f, global OSI is %2.2f', OSI,DSI,gOSI));
end

end
