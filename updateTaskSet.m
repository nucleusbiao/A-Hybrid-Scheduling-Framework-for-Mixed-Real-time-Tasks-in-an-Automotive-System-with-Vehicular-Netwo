function [taskSet, task, jobDrop] = updateTaskSet(taskSet,task, curTime, nextTime, jobProcessIndex, aetPercent)
% taskSet{i}{1} period
% taskSet{i}{2} wcetLO
% taskSet{i}{3} always 0
% taskSet{i}{4} deadline
% taskSet{i}{5} deadline

% taskSet{i}{6} aet of current job
% taskSet{i}{7} release time of current job
% taskSet{i}{8} absolute deadline of current job


jobDrop = 0;
for i = 1:size(taskSet,2)
    if isempty(taskSet{i})
        taskNum = i - 1;
        break
    end
end
if jobProcessIndex ~= -1
    for i = 1:taskNum
        if i ~= jobProcessIndex
            if taskSet{i}{8} == 0
                if ceil(curTime/taskSet{i}{1})*taskSet{i}{1} <= nextTime
                    taskSet{i}{6} = min(taskSet{i}{2}, ...
                        floor(aetPercent*taskSet{i}{2})+randi(max(1,floor((1-aetPercent)*taskSet{i}{2}))));
                    taskSet{i}{9} = -1;
                    taskSet{i}{7} = ceil(curTime/taskSet{i}{1})*taskSet{i}{1};
                    taskSet{i}{8} = ceil(curTime/taskSet{i}{1})*taskSet{i}{1} + taskSet{i}{4};

                end
            end
        else
            taskSet{i}{6} = taskSet{i}{6} - (nextTime - curTime);
            if taskSet{i}{6} ~= 0
                taskSet{i}{3} = taskSet{i}{3} + (nextTime - curTime);
            else
                taskSet{i}{8} = 0; 
                taskSet{i}{3} = 0;
            end
        end
    end
else
    for i = 1:size(taskSet, 2)
        if ~isempty(taskSet{i})
        if nextTime == ceil(curTime/taskSet{i}{1})*taskSet{i}{1}
                taskSet{i}{6} = ceil(aetPercent*taskSet{i}{2})+ceil(ceil((1-aetPercent)*taskSet{i}{2}));
%                 taskSet{i}{9} = -1;
            taskSet{i}{7} = ceil(curTime/taskSet{i}{1})*taskSet{i}{1};
            taskSet{i}{8} = ceil(curTime/taskSet{i}{1})*taskSet{i}{1} + taskSet{i}{4};
        end
        end
    end
end
for i = 1:taskNum
    if ~isempty(taskSet{i})
        if taskSet{i}{6} == 0&&taskSet{i}{9} == 1
            taskSet(i:taskNum) = taskSet(i+1:taskNum+1);
        end
    end
end
end