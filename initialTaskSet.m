function taskSet = initialTaskSet(taskSetOld, aetPercent)
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
for taskIndex = 1:size(taskSet, 2)
    if ~isempty(taskSet{taskIndex})
            taskSet{taskIndex}{6} = min(taskSet{taskIndex}{2}, ...
                floor(aetPercent*taskSet{taskIndex}{2})+randi(max(1,floor((1-aetPercent)*taskSet{taskIndex}{2}))));
            taskSet{taskIndex}{9} = -1;
        taskSet{taskIndex}{7} = 0;
        taskSet{taskIndex}{8} = taskSet{taskIndex}{4};
    end
end