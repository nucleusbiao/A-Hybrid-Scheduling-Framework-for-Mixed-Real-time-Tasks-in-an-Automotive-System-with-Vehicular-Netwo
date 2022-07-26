function jobProcessIndex = jobSearchFP(taskSetOld)
% taskSet{i}{1} period
% taskSet{i}{2} wcetLO
% taskSet{i}{3} always 0
% taskSet{i}{4} deadline
% taskSet{i}{5} deadline

% taskSet{i}{6} aet of current job
% taskSet{i}{7} release time of current job
% taskSet{i}{8} absolute deadline of current job
% taskSet{i}{9} always -1


jobProcessIndex = -1;
for i = 1:size(taskSetOld, 2)
    if ~isempty(taskSetOld{i})
    if taskSetOld{i}{8} ~= 0 
        jobProcessIndex = i;
        break       
    end
    end
end

