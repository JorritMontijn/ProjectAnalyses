vecFast = zeros(1,size(cellBehavRT,2));
vecSlow = zeros(1,size(cellBehavRT,2));
for intAnimal=1:size(cellBehavRT,2)
	vecAnim=sort(cell2mat(cellBehavRT(:,intAnimal)),'ascend');
	dblFast=mean(vecAnim(1:floor(length(vecAnim)/2)));
	dblSlow=mean(vecAnim(ceil(length(vecAnim)/2):end));
	vecFast(intAnimal) = dblFast;
	vecSlow(intAnimal) = dblSlow;
	fprintf('Anim %d; fast = %.3f;slow = %.3f\n',intAnimal,dblFast,dblSlow);
end
fprintf('Means; fast = %.3f;slow = %.3f\n',mean(vecFast),mean(vecSlow));