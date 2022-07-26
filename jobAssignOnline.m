clear

% % Load data
% -- load('folder path\data\100Task30processorCount1.1U.mat')   50 represents the number of tasks
% or you can use generate_task.m to generate new data

load('your folder path\data\100Task30processorCount1.1U.mat') 
processorCount = 30;
hardTaskCount = 50;
softTaskCount = 50;


%%  hard task assign
left = 1;
right = hardTaskCount;
middle = right;
middleTemp = 0;

% binary search   test utilization
while abs(right-left) >= 2                      
    middle0 = middle;
    taskSetTemp = hardTaskSet(1:middle, :);
    lambdaTemp = lambdaHard(:, 1:middle);
    taskCount = middle;
    numX = taskCount*processorCount;     % The number of variables
    
    f = zeros(1, numX);
    intcon = 1:numX;       %Set all variables to integers 0 and 1
    lb = zeros(numX,1);
    ub = ones(numX,1);
    
    flag2 = 1;     %flag=0, hard task can't schedule
    A = [];
    Aeq = zeros(taskCount,numX);
    for i = 1:size(lambdaTemp,1)
        A(i, ((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
        f(((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
    end
    
    for j = 1:taskCount
        for i = 1:processorCount
            Aeq(j, (i-1)*taskCount+j) = 1;
        end
    end
    b = 0.693*ones(processorCount,1);
    beq = ones(taskCount,1);
    options = optimoptions('intlinprog','Display','off');
    
    [xTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);  % Solve the 01 programming problem
    if ~isempty(xTemp)
        x = xTemp;
        middleTemp = middle;
        left = middle;
    else
        right = middle;
    end
    middle = ceil((left+right)/2);
end
middle = middle0;
if middle == hardTaskCount
    disp('success for hard real-time tasks');
    flag1 = 1;
else
    flag1 = 0;
end
for i = 1:processorCount
    ecuTaskSet{i} = [];
end
if ~isempty(xTemp)
    x = xTemp;
    for i = 1:processorCount
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2
                ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
            end
        end
    end
else
    middle = middleTemp;
    taskSetTemp = hardTaskSet(1:middle, :);
    lambdaTemp = lambdaHard(:, 1:middle);
    taskCount = middle;
    numX = taskCount*processorCount;
    for i = 1:processorCount
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2
                ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
            end
        end
    end
end
ecuUtilizationTemp = [];
for i = 1:processorCount
    if ~isempty(ecuTaskSet{i})
        ecuUtilizationTemp(i) = sum(ecuTaskSet{i}(:,4) .*ecuTaskSet{i}(:,3)./ecuTaskSet{i}(:,1));  %Calculate the utilization used by each processor
    else
        ecuUtilizationTemp(i) = 0;
    end
end

%% When utilization tests cannot allocate all hard tasks, supply and demand analysis is used to allocate the remaining hard tasks. 
%  If all hard tasks are assigned through Demand-supply analysis, then soft tasks are directly assigned by Demand-supply analysis.

% Assign soft tasks when hard tasks have already been assigned
if flag1 == 0
    hardTaskSetLeft = hardTaskSet(middle+1:end, :);
    lambdaHardLeft = lambdaHard(:, middle+1:end);
    UnassignedTasks = [];
    for j = 1:size(hardTaskSetLeft, 1)

        % Sort processors by utilization, with the largest remaining utilization taking precedence
        [result, index] = sort(ecuUtilizationTemp, 'ascend'); 
        for i = 1:processorCount
            processorID = index(i);
            taskCountThisProcessor = size(ecuTaskSet{processorID}, 1);

            % Perform unsimplified utilization tests
            if ecuUtilizationTemp(processorID) + lambdaHardLeft(processorID,j)*hardTaskSetLeft(j,3)/hardTaskSetLeft(j,1) < (taskCountThisProcessor+1)*(2^(1/(taskCountThisProcessor+1))-1)
                ecuUtilizationTemp(processorID) = ecuUtilizationTemp(processorID) + lambdaHardLeft(processorID,j)*hardTaskSetLeft(j,3)/hardTaskSetLeft(j,1);
                ecuTaskSet{processorID} = [ecuTaskSet{processorID}; hardTaskSetLeft(j,:), lambdaHardLeft(processorID,j)];
                display('success for this task by util');
                break
            else
                % Demand-supply analysis
                taskArray = [ecuTaskSet{processorID};hardTaskSetLeft(j,:), lambdaHardLeft(processorID,j)];
                schFlag = schTestSupplyDemand(taskArray);
                if schFlag
                    ecuUtilizationTemp(processorID) = ecuUtilizationTemp(processorID) + lambdaHardLeft(processorID,j)*hardTaskSetLeft(j,3)/hardTaskSetLeft(j,1);
                    ecuTaskSet{processorID} = [ecuTaskSet{processorID}; hardTaskSetLeft(j,:), lambdaHardLeft(processorID,j)];
                    display('success for this task by mpa');
                    break
                end
            end
            
            if eq(i, processorCount)
                display('failure for this task');
                UnassignedTasks = [UnassignedTasks; hardTaskSetLeft(j,:), lambdaHardLeft(:,j)'];
            end
        end
    end
    if ~isempty(UnassignedTasks)
        flag2 = 0;
    end
end

%% assign soft task 
if flag1 == 1
    left = 1;
    right = softTaskCount;
    middleTemp = 0;
    while abs(right-left) >= 2
        middle = ceil((left+right)/2);
        taskSetTemp = softTaskSet(1:middle, :);
        lambdaTemp = lambdaSoft(:, 1:middle);
        taskCount = middle;
        numX = taskCount*processorCount;
        Aeq = zeros(taskCount,numX);
        f = zeros(1, numX);
        intcon = 1:numX;
        lb = zeros(numX,1);
        ub = ones(numX,1);
        
        A = [];
        Aeq = zeros(taskCount,numX);
        for i = 1:size(lambdaTemp,1)
            A(i, ((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
            f(((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
        end
        for j = 1:taskCount
            for i = 1:processorCount
                Aeq(j, (i-1)*taskCount+j) = 1;
            end
        end
        b = 0.693*ones(processorCount,1) - ecuUtilizationTemp';
        beq = ones(taskCount,1);
        [yTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);
        if ~isempty(yTemp)
            y = yTemp;
            middleTemp = middle;
            left = middle;
        else
            right = middle;
        end
    end
    if ~isempty(yTemp)
        for i = 1:processorCount
            for j = 1:taskCount
                if abs(y((i-1)*taskCount+j)-1) < 1e-2
                    ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
                end
            end
        end
    else
        middle = middleTemp;
        taskSetTemp = softTaskSet(1:middle, :);
        lambdaTemp = lambdaSoft(:, 1:middle);
        taskCount = middle;
        numX = taskCount*processorCount;
        for i = 1:processorCount
            for j = 1:taskCount
                if abs(y((i-1)*taskCount+j)-1) < 1e-2
                    ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
                end
            end
        end
    end
    softTaskSetLeft = [];
    if middle+1~=softTaskCount
        softTaskSetLeft = softTaskSet(middle+1:end, :);
    end
    ecuUtilization = [];
    for i = 1:processorCount
        if ~isempty(ecuTaskSet{i})
            ecuUtilization(i) = sum(ecuTaskSet{i}(:,4) .*ecuTaskSet{i}(:,3)./ecuTaskSet{i}(:,1));
        end
    end
    lambdaSoftLeft = lambdaSoft(:, middle+1:end);
else
    if flag2 == 1
        lambdaSoftLeft = lambdaSoft;
        softTaskSetLeft = softTaskSet;
        ecuUtilization = ecuUtilizationTemp;
    end
end

%% If flag1 or flag2 is 1, the hard task is fully scheduled; 
% flag2=1, soft tasks are directly assigned by Demand-supply analysis.

if flag1 == 1||flag2 == 1
    UnassignedTasks = [];
    for j = 1:size(softTaskSetLeft, 1)
        [result, index] = sort(ecuUtilization, 'ascend');
        for i = 1:processorCount
            processorID = index(i);
            taskCountThisProcessor = size(ecuTaskSet{processorID}, 1);
            if ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1) < (taskCountThisProcessor+1)*(2^(1/(taskCountThisProcessor+1))-1)
                ecuUtilization(processorID) = ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1);
                ecuTaskSet{processorID} = [ecuTaskSet{processorID}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
                display('success for this task by util');
                break
            else
                taskArray = [ecuTaskSet{processorID}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
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
        end
    end
end

%% assign online
taskSet = {}; 
for n = 1:processorCount
    ecuTaskSet{n} = sortrows(ecuTaskSet{n}, 1);
end
for n = 1:processorCount
    for i = 1:25
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
        task1{1,1} = ecuTaskSet{n}(i,1);  
        task1{2,1} = ecuTaskSet{n}(i,3);  
        task1{3,1} = 0;
        task1{4,1} = ecuTaskSet{n}(i,2);
        task1{5,1} = ecuTaskSet{n}(i,4);
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
timeLength = 10^6;
curTime = 0;
aetPercent =  0.6;
% 
for i = 1:size(UnassignedTasks,1)
    task{1,i} = UnassignedTasks(i,1);  %
    task{2,i} = UnassignedTasks(i,3);  %
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
lambdaLeft = UnassignedTasks(:,4:end);
%%
Ect = zeros(1000,1);
ns = 0;nsc = 0;nse = 0;
nf = 0;nt = 0;
while curTime <= timeLength
    num =1;
    for i = 1:processorCount
        [nextTimes(i), jobProcessIndex(i)] = CalculateNextTime(taskSet(i,:), ...
            curTime);
    end
    taskTemp = sortrows(task',7)';            %
    for i = 2:size(task,2)
        if taskTemp{7,i} ~= taskTemp{7,1}
            break
        end
        num = num + 1;           %
    end
    nextTime = min([nextTimes taskTemp{7,1}]);
    if nextTime == taskTemp{7,1}
        for i = 1:num
            t1 = clock;
            [schFlag,processorID,taskSet,pos] = schTestOnline(taskSet,taskTemp(:,i),lambdaLeft,taskTemp{7,i});
            if schFlag
                ns = ns + 1;
                disp('success for this job');
                disp(['taskID = ',num2str(taskTemp{10,i}),' ','processorID = ',num2str(processorID)]);
            else
                if randi([0 1])
                    if rand > 0.5
%                         disp('job assigned to Cloud');
                        nsc = nsc + 1;
                    else
%                         disp('failure for this job');
                        nf = nf + 1;
                    end
                else
                    if rand > 0.5
%                         disp('job assigned to Edge');
                        nse = nse + 1;
                    else
%                         disp('failure for this job');
                        nf = nf + 1;
                    end
                end
            end
            t2 = clock;
            nt = nt + 1;
            Ect(nt) = etime(t2,t1);
            task{7,taskTemp{10,i}} = task{7,taskTemp{10,i}} + task{1,taskTemp{10,i}};
            task{8,taskTemp{10,i}} = task{7,taskTemp{10,i}} + task{1,taskTemp{10,i}};
        end
    end
    if nt>999
        break
    end
    if nextTime ~= curTime
        for i = 1:processorCount
            [taskSet(i,:), task, jobDrop] = updateTaskSet(taskSet(i,:), task, curTime, nextTime, jobProcessIndex(i), aetPercent); %更新
        end
    end
    curTime = nextTime;
end

