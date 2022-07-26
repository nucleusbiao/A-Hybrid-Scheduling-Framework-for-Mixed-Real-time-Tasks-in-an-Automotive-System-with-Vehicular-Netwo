%  This code is the same as USDPLUSE.m
% With less data, you can run directly to see the results

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
    hardTaskSet(i,:) = [period, period, 200+10*randi(10)];
end

for i = 1:softTaskCount
    period = periodArray(randi(9));
    softTaskSet(i,:) = [period, period, 200+10*randi(10)];
end

lambdaHard = (5+randi(10, processorCount, hardTaskCount))/10;
lambdaSoft = (5+randi(10, processorCount, softTaskCount))/10;

hardTaskSet = sortrows(hardTaskSet, 1);
softTaskSet = sortrows(softTaskSet, 1);
%%
numX = hardTaskCount*processorCount;
f = ones(1, numX);
intcon = 1:numX;
lb = zeros(numX,1);
ub = ones(numX,1);
%%
exitflag = 1;
left = 1;
right = softTaskCount;

while abs(right-left) > 2
    middle = ceil((left+right)/2);
    softTaskSetTemp = softTaskSet(1:middle, :);
    lambdaSoftTemp = lambdaSoft(:, 1:middle);
    
    taskSetTemp = [hardTaskSet; softTaskSetTemp];
    lambdaTemp = [lambdaHard, lambdaSoftTemp];
    
    taskCount = hardTaskCount + middle;
    numX = taskCount*processorCount;
    f = ones(1, numX);
    intcon = 1:numX;
    lb = zeros(numX,1);
    ub = ones(numX,1);

    A = [];
    Aeq = zeros(taskCount, numX);
    for i = 1:size(lambdaTemp,1)
        A(i, ((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
        % f(((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
    end

    for j = 1:taskCount
        for i = 1:processorCount
            Aeq(j, (i-1)*taskCount+j) = 1;
        end
    end

    b = 0.693*ones(processorCount,1);
    beq = ones(taskCount,1);
    options = optimoptions('intlinprog','Display','off');
    [xTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);
    if ~isempty(xTemp)
        x = xTemp;
        middleTemp = middle;
        left = middle;
    else
        right = middle;
    end        
end

for i = 1:processorCount
    ecuTaskSet{i} = [];
end
    
if ~isempty(xTemp)
    for i = 1:processorCount
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2
                ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
            end
        end
    end
else
    middle = middleTemp;
    softTaskSetTemp = softTaskSet(1:middle, :);
    lambdaSoftTemp = lambdaSoft(:, 1:middle);
    
    taskSetTemp = [hardTaskSet; softTaskSetTemp];
    lambdaTemp = [lambdaHard, lambdaSoftTemp];
    taskCount = hardTaskCount + middle;
    
    for i = 1:processorCount
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2
                ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
            end
        end
    end
end

softTaskSetLeft = softTaskSet(middle+1:end, :);
lambdaSoftLeft = lambdaSoft(:, middle+1:end);
%%
ecuUtilization = [];
for i = 1:processorCount
    ecuUtilization(i) = sum(ecuTaskSet{i}(:,4) .*ecuTaskSet{i}(:,3)./ecuTaskSet{i}(:,1));
end
%%
UnassignedTasks = [];
for j = 1:size(softTaskSetLeft, 1)
    [result, index] = sort(ecuUtilization, 'ascend');
    for i = 1:processorCount
        processorID = index(i);
        taskCountThisProcessor = size(ecuTaskSet{i}, 1);
        if ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1) < (taskCountThisProcessor+1)*(2^(1/(taskCountThisProcessor+1))-1)
            ecuUtilization(processorID) = ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1);
            ecuTaskSet{i} = [ecuTaskSet{i}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
            display('success for this task by util');
            break
        else
            taskArray = [ecuTaskSet{i}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
            schFlag = schTestSupplyDemand(taskArray);
            if schFlag
                ecuUtilization(processorID) = ecuUtilization(processorID) + lambdaSoftLeft(processorID,j)*softTaskSetLeft(j,3)/softTaskSetLeft(j,1);
                ecuTaskSet{i} = [ecuTaskSet{i}; softTaskSetLeft(j,:), lambdaSoftLeft(processorID,j)];
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


