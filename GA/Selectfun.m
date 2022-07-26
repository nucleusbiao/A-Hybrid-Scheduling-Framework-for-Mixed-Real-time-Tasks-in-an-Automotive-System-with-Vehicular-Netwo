function [newgroup] = Selectfun(group,fit,processorCount,hardTaskCount)
fitsum = sum(fit);
newgroup = zeros(size(group,1),size(group,2));
[maxfit,index] = max(fit);
for i = 1:size(group,2)
    p(i) = fit(i)/fitsum;
end
p(index) = 0;

for i = 1:size(group,2)
    r = rand();
    q = 0;    
    for j = 1:size(fit,2)
        q = q + p(j);
        if q >= r
            newgroup(:,i) = group(:,j);
            p(j) = 0;
            break
        end       
        newgroup(1:hardTaskCount,i) = processorCount/10*rand(hardTaskCount,1);             
        newgroup(hardTaskCount+1:end,i) = (processorCount + 15)/10*rand(size(group,1)-hardTaskCount,1);        
    end
    if i == index
        newgroup(:,i) = group(:,index);
    end
end
end

