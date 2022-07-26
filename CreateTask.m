function [hardTaskSet,softTaskSet,lambdaHard,lambdaSoft] = CreateTask(processorCount,hardTaskCount,softTaskCount)
% generate task
%{
    If you want to generate multiple sets of data to use in 
    USD.m USDPLUS.m EXACT.m,  use generate_task.mlx
%}

periodArray1 = [20 40 80 200 400 800];   % hard task period
periodArray2 = [25 50 100 250 500 1000];            % soft task period
taskCount = hardTaskCount + softTaskCount;

while 1
    vectU = UUniFast(taskCount,16);          % Task utilization 
    if max(vectU)<1
        break
    end
end
for m = 1:hardTaskCount
    period = periodArray1(randi(6));
    hardTaskSet(m,:) = [period, period, ceil(vectU(m)*period)];      % period  deadline  worst execution time
end

for n = 1:softTaskCount
    period = periodArray2(randi(6));
    softTaskSet(n,:) = [period, period, ceil(vectU(m+hardTaskCount)*period)];
end

lambdaHard = (8+randi(4, processorCount, hardTaskCount))/10;             % The execution efficiency of the processor
lambdaSoft = (8+randi(4, processorCount, softTaskCount))/10;
hardTaskSet = sortrows(hardTaskSet, 1);
softTaskSet = sortrows(softTaskSet, 1);

end

