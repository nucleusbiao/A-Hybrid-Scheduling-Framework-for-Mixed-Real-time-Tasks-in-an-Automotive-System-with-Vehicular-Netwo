% % clear

hardTaskSet0 = [];
softTaskSet0 = [];


% % Load data
% -- load('folder path\data\50.mat')   50 represents the number of tasks
% or you can use generate_task.m to generate new data

load('your folder path\data\50.mat') 
for k = 1:10

    hardTaskSet0 = hardTaskSet{k};
    softTaskSet0 = softTaskSet{k};
    lambdaHard0 = lambdaHard{k}(1:processorCount,:);
    lambdaSoft0 = lambdaSoft{k}(1:processorCount,:);
    
    hardTaskSet0 = sortrows(hardTaskSet0, 1);
    softTaskSet0 = sortrows(softTaskSet0, 1);
    %% Hard Real-Time Task Assignment

    lambdaTemp = lambdaHard0;
    taskSetTemp = hardTaskSet0;
    taskCount = hardTaskCount;
    numX = taskCount*processorCount;
    f = zeros(1, numX);
    intcon = 1:numX;
    lb = zeros(numX,1);
    ub = ones(numX,1);
    flag = 1;
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
    t1 = clock;
    [xTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);
    
    for i = 1:processorCount
        ecuTaskSet{i} = [];
    end

    % Decode and generate the result
    if ~isempty(xTemp)
        disp('success for hard real-time tasks');
        x = xTemp;
        for i = 1:processorCount
            for j = 1:taskCount
                if abs(x((i-1)*taskCount+j)-1) < 1e-2
                    ecuTaskSet{i} = [ecuTaskSet{i}; taskSetTemp(j,:), lambdaTemp(i,j)];
                end
            end
        end
    else
        disp('failure');
        flag = 0;
    end
    ecuUtilizationTemp = [];
    for i = 1:processorCount
        if ~isempty(ecuTaskSet{i})
            ecuUtilizationTemp(i) = sum(ecuTaskSet{i}(:,4) .*ecuTaskSet{i}(:,3)./ecuTaskSet{i}(:,1));
        end
    end
    %% binary search to assign Soft Real-Time Task
    if flag == 1
        left = 1;
        right = softTaskCount;
        middleTemp = 0;
        while abs(right-left) >= 2
            middle = ceil((left+right)/2);
            taskSetTemp = softTaskSet0(1:middle, :);
            lambdaTemp = lambdaSoft0(:, 1:middle);
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
            taskSetTemp = softTaskSet0(1:middle, :);
            lambdaTemp = lambdaSoft0(:, 1:middle);
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
            softTaskSetLeft = softTaskSet0(middle+1:end, :);
        end
        ecuUtilization = [];
        for i = 1:processorCount
            if ~isempty(ecuTaskSet{i})
                ecuUtilization(i) = sum(ecuTaskSet{i}(:,4) .*ecuTaskSet{i}(:,3)./ecuTaskSet{i}(:,1));
            end
        end
        lambdaSoftLeft = lambdaSoft0(:, middle+1:end);
        
        
        %% Demand-supply schedulability test to assign the rest of the tasks
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
    if flag == 1
        USDtaskN(k) = hardTaskCount + softTaskCount - size(UnassignedTasks,1);
    else
        USDtaskN(k) = 0;
    end
    t2 = clock;
    USDruntime(k) = etime(t2,t1);
end
% save('C:\Users\BUCT\Desktop\新建文件夹\1\数据\USDruntime','USDruntime')
% save('C:\Users\BUCT\Desktop\新建文件夹\1\数据\USDtaskN','USDtaskN')
% disp('fininsh!!!!!!!!!!')
