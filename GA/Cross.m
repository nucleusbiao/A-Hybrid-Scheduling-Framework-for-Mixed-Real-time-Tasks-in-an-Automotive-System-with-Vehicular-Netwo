function [group] = Cross(group,A,processorCount)
fit = fitness(group,A,processorCount);
[maxfit,index] = max(fit);
bestsingle = group(:,index);
for i = 1:2:size(group,2)
    pos = ceil(rand()*size(group,1));
    
    groupTemp = group(pos:size(group,1),i);
    
    group(pos:size(group,1),i) = group(pos:size(group,1),mod(i,size(group,2))+1);
    group(pos:size(group,1),mod(i,size(group,2))+1) = groupTemp;
end
group(:,1) = bestsingle;
end

