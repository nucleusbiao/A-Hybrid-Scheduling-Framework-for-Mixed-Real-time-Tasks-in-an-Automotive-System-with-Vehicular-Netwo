function [nextTime, jobProcessIndex] = CnextTime(taskSetOld, ...
    curTime)
% taskSet{i}{1} period
% taskSet{i}{2} wcetLO
% taskSet{i}{3} always 0
% taskSet{i}{4} deadline
% taskSet{i}{5} deadline

% taskSet{i}{6} aet of current job
% taskSet{i}{7} release time of current job
% taskSet{i}{8} absolute deadline of current job
% taskSet{i}{9} always -1

taskSet = taskSetOld;
jobProcessIndex = jobSearchFP(taskSetOld);
if jobProcessIndex ~= -1
        nextTime = curTime + taskSet{jobProcessIndex}{6};
    for i = 1:jobProcessIndex-1
        if ceil(curTime/taskSet{i}{1})*taskSet{i}{1} < nextTime
            nextTime = ceil(curTime/taskSet{i}{1})*taskSet{i}{1};
        end
    end
else
    nextTime = 1e12;
    for i = 1:size(taskSet, 2)
        if ~isempty(taskSet{i})
        if taskSet{i}{8} == 0
            nextTimeTemp = (floor(curTime/taskSet{i}{1})+1)*taskSet{i}{1};
            if nextTimeTemp < nextTime
                nextTime = nextTimeTemp;
            end
        end
        end
    end
end
