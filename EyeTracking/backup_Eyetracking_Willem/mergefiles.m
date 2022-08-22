%Find all video part files and load them in 1 variable loadedData

fileList = dir([filename '_video_*']);
nFiles = length(fileList);

%# loop through files and write results into loadedData
for iFile = 1:nFiles
    tmp = load(fileList(iFile).name);
    %loadedData(iFile,:) = struct2cell(tmp);
    loadedData(iFile) = tmp;
end

%Merge structs in loaddata

fields = fieldnames(loadedData(1).video_parts)';

% %Nodig voor maken van cellfun
% for i = 1:nFiles
%     
%     if i == 1
%         stringloadedData = ['[loadedData(1,' num2str(i) ').video_parts.(f) '] ;
%     elseif i == nFiles
%         stringloadedData = [stringloadedData 'loadedData(1,' num2str(i) ').video_parts.(f)]'] ;
%     else
%         stringloadedData = [stringloadedData 'loadedData(1,' num2str(i) ').video_parts.(f) '] ;
%     end
%     
% end

B = cat(3,loadedData.video_parts);
out = reshape(B(1,1,:),1,[]);
fields(2,:) = cellfun(@(f) [out(1,:).(f)], fields, 'unif', false);

%fields(2,:) = cellfun(@(f) [loadedData(1,1).video_parts.(f) loadedData(1,2).video_parts.(f) loadedData(1,3).video_parts.(f) loadedData(1,4).video_parts.(f) loadedData(1,5).video_parts.(f) loadedData(1,6).video_parts.(f) loadedData(1,7).video_parts.(f) loadedData(1,8).video_parts.(f) loadedData(1,9).video_parts.(f) loadedData(1,10).video_parts.(f) loadedData(1,11).video_parts.(f) loadedData(1,12).video_parts.(f) loadedData(1,13).video_parts.(f) loadedData(1,14).video_parts.(f) loadedData(1,15).video_parts.(f) loadedData(1,16).video_parts.(f) loadedData(1,17).video_parts.(f) loadedData(1,18).video_parts.(f) loadedData(1,19).video_parts.(f) loadedData(1,20).video_parts.(f) loadedData(1,21).video_parts.(f) loadedData(1,22).video_parts.(f) loadedData(1,23).video_parts.(f) loadedData(1,24).video_parts.(f) loadedData(1,25).video_parts.(f) loadedData(1,26).video_parts.(f) loadedData(1,27).video_parts.(f)], fields, 'unif', false);
%fields(2,:) = cellfun(@(f) cat(3,loadedData.(f).video_parts), fields, 'unif', false);

video_heel = struct(fields{:});
