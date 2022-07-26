clear
periodArray = [5000 800 1000 1600 500 600 900 1000 1000];
wcetLOArray = [20 20 20 20 20 200 200 200 200];

processorCount = 10;


hardTaskSet = [];
softTaskSet = [];
hardTaskCount = 25;
softTaskCount = 30;


for i = 1:hardTaskCount
    period = periodArray(randi(9));                                    
    hardTaskSet(i,:) = [period, period, 200+10*randi(10)];      % period  deadline  worst execution time  
end

for i = 1:softTaskCount
    period = periodArray(randi(9));
    softTaskSet(i,:) = [period, period, 150+10*randi(10)];
end

lambdaHard = (5+randi(10, processorCount, hardTaskCount))/10;   % The execution efficiency of the processor
lambdaSoft = (5+randi(10, processorCount, softTaskCount))/10;

hardTaskSet = sortrows(hardTaskSet, 1);
softTaskSet = sortrows(softTaskSet, 1);

%%
exitflag = 1;
    
    taskSetTemp = [hardTaskSet; softTaskSet];
    lambdaTemp = [lambdaHard, lambdaSoft];
    
    taskCount = hardTaskCount + softTaskCount;
    numX = taskCount*(processorCount+1);
    f = zeros(1, numX);
    intcon = 1:numX;
    lb = zeros(numX,1);
    ub = ones(numX,1);
    
    for j=1:processorCount
        f((j-1)*taskCount+1:j*taskCount) = -1;  
    end
    A = zeros(processorCount, numX);
    Aeq = zeros(taskCount, numX);
    for i = 1:size(lambdaTemp,1)
        A(i, ((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
        % f(((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
    end

    for j = 1:taskCount
        if j < hardTaskCount + 1
            for i = 1:processorCount
                Aeq(j, (i-1)*taskCount+j) = 1;
            end
        else
            for i = 1:processorCount+1
                Aeq(j, (i-1)*taskCount+j) = 1;
            end
        end
    end

    b = 0.693*ones(processorCount,1);
    beq = ones(taskCount,1);
    options = optimoptions('intlinprog','Display','off');
    [xTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);
    x=xTemp;

for i = 1:processorCount
    ecuTaskSet{i} = [];
end
 virtprocessor=[];    % Virtual processor
if ~isempty(xTemp)
    for i = 1:processorCount+1
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2&&i<processorCount+1
                ecuTaskSet{i} = [ecuTaskSet{i}; j, taskSetTemp(j,:), lambdaTemp(i,j)];
            end
            if abs(x((i-1)*taskCount+j)-1) < 1e-2&&i==processorCount+1
                virtprocessor=[virtprocessor; j, taskSetTemp(j,:)];
            end
        end
    end
end
softTaskSetLeft = virtprocessor;
lambdaSoftLeft = zeros(processorCount,size(virtprocessor,1));
for i = 1:size(virtprocessor,1)
    lambdaSoftLeft(:,i) = lambdaTemp(:,virtprocessor(i,1));
end
ecuUtilization = [];
for i = 1:processorCount
    ecuUtilization(i) = sum(ecuTaskSet{i}(:,5) .*ecuTaskSet{i}(:,4)./ecuTaskSet{i}(:,2));
end

%% 

%% Unassigned soft task
UnassignedTasks = [];
for j = 1:size(softTaskSetLeft, 1)
    [result, index] = sort(ecuUtilization, 'ascend');
    for i = 1:processorCount
        processorID = index(i);
        taskCountThisProcessor = size(ecuTaskSet{processorID}, 1);
        if ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1) < (taskCountThisProcessor+1)*(2^(1/(taskCountThisProcessor+1))-1)
            ecuUtilization(processorID) = ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1);
            ecuTaskSet{processorID} = [ecuTaskSet{processorID};softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
            display('success for this task by util');
            break
        else
            taskArray = [ecuTaskSet{processorID}(:,2:end);softTaskSetLeft(j,2:end), lambdaSoftLeft(processorID,j)];%softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)
            schFlag = schTestSupplyDemand(taskArray);
            if schFlag
                ecuUtilization(processorID) = ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1);
                ecuTaskSet{processorID} = [ecuTaskSet{processorID}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
                display('success for this task by mpa');
                break
            end
        end
        
        if eq(i, processorCount)
            display('failure for this task');
            UnassignedTasks = [UnassignedTasks; softTaskSetLeft(j,:), lambdaSoftLeft(:,j)'];
        end
% for i = 1:size(softTaskSetLeft,1)
%      UnassignedTasks = [UnassignedTasks; softTaskSetLeft(j,:), lambdaSoftLeft(:,j)'];
% end
    end
end
%% assign online
taskSet = {};
for n = 1:processorCount
    for i = 1:10
        taskSet{n,i}={};
    end
end
for n = 1:processorCount
% taskSet{i}{1} period
% taskSet{i}{2} wcetLO
% taskSet{i}{3} always 0
% taskSet{i}{4} deadline
% taskSet{i}{5} deadline

% taskSet{i}{6} aet of current job
% taskSet{i}{7} release time of current job
% taskSet{i}{8} absolute deadline of current job
% taskSet{i}{9} 
    for i = 1:size(ecuTaskSet{n}, 1)
        task1{1,1} = ecuTaskSet{n}(i,2);  % period
        task1{2,1} = ecuTaskSet{n}(i,4);  % worst execution time
        task1{3,1} = 0;
        task1{4,1} = ecuTaskSet{n}(i,3);
        task1{5,1} = ecuTaskSet{n}(i,5);
        task1{6,1} = 0;
        task1{7,1} = 0;
        task1{8,1} = 0;
        task1{9,1} = 0;
        task1{10,1} = 0;
        taskSet{n,i} = task1;
    end
end
lambdainsertTask = [];

%%
timeLength = 10000;
curTime = 0;
aetPercent = 0.2;
for i = 1:size(UnassignedTasks,1)
    task{1,i} = UnassignedTasks(i,2); 
    task{2,i} = UnassignedTasks(i,4);  
    task{3,i} = 0;
    task{4,i} = UnassignedTasks(i,2);
    task{5,i} = 1;
    task{6,i} = min(task{2,i}, ...
        floor(aetPercent*task{2,i})+randi(max(1,floor((1-aetPercent)*task{2,i}))));
    task{7,i} = 0;
    task{8,i} = task{7,i} + task{4,i};
    task{9,i} = 1;
    task{10,i} = i;
end
for i = 1:processorCount
    taskSet(i,:) = initialTaskSet(taskSet(i,:), aetPercent);
end
lambdaSoftLeft = UnassignedTasks(:,5:end);
%%
clc
n = 0;
n0 = 0;
while curTime <= timeLength
    num =1;
    for i = 1:processorCount
        [nextTimes(i), jobProcessIndex(i)] = CnextTime(taskSet(i,:), ...
            curTime);
    end
    taskTemp = sortrows(task',7)';
    for i = 2:size(task,2)
        if taskTemp{7,i} ~= taskTemp{7,1}
            break
        end
        num = num + 1;
    end
    nextTime = min([nextTimes taskTemp{7,1}]);
    if nextTime == taskTemp{7,1}
        for i = 1:num
            [schFlag,processorID,taskSet,pos] = schTestOnline(taskSet,taskTemp(:,i),lambdaSoftLeft,taskTemp{7,i});
            if schFlag
                disp('success for this job');
                disp(['taskID = ',num2str(taskTemp{10,i}),' ','processorID = ',num2str(processorID)]);
                disp(pos)
                n0 = n0 + 1;
            else
%                 disp('failure for this task');
                n = n + 1;
            end
            task{7,taskTemp{10,i}} = task{7,taskTemp{10,i}} + task{1,taskTemp{10,i}};
            task{8,taskTemp{10,i}} = task{7,taskTemp{10,i}} + task{1,taskTemp{10,i}};
        end
    end
    if nextTime ~= curTime
        for i = 1:processorCount
            [taskSet(i,:), task, jobDrop] = updateTaskSet(taskSet(i,:), task, curTime, nextTime, jobProcessIndex(i), aetPercent);
        end
    end
    curTime = nextTime;
end



