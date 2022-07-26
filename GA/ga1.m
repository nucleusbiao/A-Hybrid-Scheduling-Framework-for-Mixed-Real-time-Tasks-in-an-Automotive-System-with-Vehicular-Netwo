function [x,fval] = ga1(hardtaskCount,softtaskCount,processorCount,population,A,MaxIterations,p)

taskCount = hardtaskCount + softtaskCount;
while 1
[group] = initgroup(population,hardtaskCount,softtaskCount,processorCount);
fit = fitness(group,A,processorCount);
if max(fit)~=0
    break
end
end

for i = 1:MaxIterations
    fit = fitness(group,A,processorCount);
    max(fit)
    if max(fit) == taskCount
        break
    end
    [group] = Selectfun(group,fit,processorCount,hardtaskCount);
    [group] = Cross(group,A,processorCount);
    fit = fitness(group,A,processorCount);
    [maxfit,index] = max(fit);
    for j = 1:size(group,2)
        if rand()<p&&j~=index
            group(:,j) = mutation(group(:,j),processorCount);
        end
    end
end
[fval,index] = max(fit);
xt = group(:,index);
x = zeros(taskCount*(processorCount+1),1);
for i = 1:size(xt,1)
    if xt(i)<processorCount/10
        x((ceil(xt(i)*10)-1)*taskCount+i) = 1;
    else
        x(processorCount*taskCount+i) = 1;
    end
end
end

