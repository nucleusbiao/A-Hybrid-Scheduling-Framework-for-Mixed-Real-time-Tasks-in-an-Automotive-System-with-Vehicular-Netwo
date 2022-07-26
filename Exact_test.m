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

lambdaHard = (5+randi(10, processorCount, hardTaskCount))/10;  % The execution efficiency of the processor
lambdaSoft = (5+randi(10, processorCount, softTaskCount))/10;

hardTaskSet = sortrows(hardTaskSet, 1);
softTaskSet = sortrows(softTaskSet, 1);

%%
exitflag = 1;
left = 1;
right = softTaskCount;

    taskSetTemp = [hardTaskSet; softTaskSet];
    lambdaTemp = [lambdaHard, lambdaSoft];
    
    taskCount = hardTaskCount + softTaskCount;
    numX = taskCount*(processorCount+1);
    f = zeros(1, numX);
    intcon = 1:numX;
    lb = zeros(numX,1);
    ub = ones(numX,1);
    
    for j=1:processorCount
        f((j-1)*taskCount+1:j*taskCount) = -1;   %lambdaTemp(j,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
    end
    A = zeros(processorCount, numX);
    Aeq = zeros(taskCount, numX);
    for i = 1:size(lambdaTemp,1)
        A(i, ((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
        % f(((i-1)*taskCount+1):i*taskCount) = lambdaTemp(i,:) .* (taskSetTemp(:,3)./taskSetTemp(:,1))';
    end

    for j = 1:taskCount
        for i = 1:processorCount+1
            Aeq(j, (i-1)*taskCount+j) = 1;
        end
    end
    b = 0.693*ones(processorCount,1);
    beq = ones(taskCount,1);
    options = optimoptions('intlinprog','Display','off');

    % Solve the 01 programming problem
    [xTemp,fval,exitflag,output] = intlinprog(f,intcon,A,b, Aeq, beq, lb, ub, options);
    x = xTemp;


for i = 1:processorCount
    ecuTaskSet{i} = [];
end
 virtprocessor=[];   
if ~isempty(xTemp)
    for i = 1:processorCount+1
        for j = 1:taskCount
            if abs(x((i-1)*taskCount+j)-1) < 1e-2&&i<processorCount+1
                ecuTaskSet{i} = [ecuTaskSet{i}; j, taskSetTemp(j,:), lambdaTemp(i,j)];
            end
            if abs(x((i-1)*taskCount+j)-1) < 1e-2&&i==processorCount+1
                virtprocessor = [virtprocessor; j, taskSetTemp(j,:)];
            end
        end
    end
end

