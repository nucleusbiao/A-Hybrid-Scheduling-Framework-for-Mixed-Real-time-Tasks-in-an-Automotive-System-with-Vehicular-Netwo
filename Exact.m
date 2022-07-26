clear
% clear


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

%% Exact approach

left = 1;
right = softTaskCount;

taskSetTemp = [hardTaskSet0; softTaskSet0];
lambdaTemp = [lambdaHard0, lambdaSoft0];

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
t1=clock;

% Solve the 01 programming problem
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
    ExacttaskN(k) = hardTaskCount + softTaskCount - size(virtprocessor,1);
else
    ExacttaskN(k) = 0;
end
softTaskSetLeft = virtprocessor;
lambdaSoftLeft = zeros(processorCount,size(virtprocessor,1));
for i = 1:size(virtprocessor,1)
    lambdaSoftLeft(:,i) = lambdaTemp(:,virtprocessor(i,1));
end
ecuUtilization = [];
for i = 1:processorCount
    if ~isempty(ecuTaskSet{i})
        ecuUtilization(i) = sum(ecuTaskSet{i}(:,5) .*ecuTaskSet{i}(:,4)./ecuTaskSet{i}(:,2));
    end
end
t2 = clock;

Exactruntime(k) = etime(t2,t1);

end
% save('C:\Users\BUCT\Desktop\新建文件夹\1\数据\Exact-runtime','Exactruntime')
% save('C:\Users\BUCT\Desktop\新建文件夹\1\数据\Exact-taskN','ExacttaskN')
% disp('finish!!!!!!!!!!')
%% 

            
