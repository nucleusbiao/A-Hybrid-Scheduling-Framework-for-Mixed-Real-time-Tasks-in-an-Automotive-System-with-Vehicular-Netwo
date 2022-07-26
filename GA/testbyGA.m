clear
periodArray = [20 25 40 50 80 100 200 250 400 500 800 1000];

processorCount = 20;
hardTaskSet = [];
softTaskSet = [];
hardTaskCount = 35;
softTaskCount = 35;
vectU = UUniFast(40,9);
for i = 1:hardTaskCount
    period = periodArray(randi(12));                                    
    hardTaskSet(i,:) = [period, period, ceil(vectU(randi(40))*period)];      %
end

for i = 1:softTaskCount
    period = periodArray(randi(12));
    softTaskSet(i,:) = [period, period, ceil(vectU(randi(40))*period)];
end

lambdaHard = (5+randi(10, processorCount, hardTaskCount))/10;
lambdaSoft = (5+randi(10, processorCount, softTaskCount))/10;

hardTaskSet = sortrows(hardTaskSet, 1);
softTaskSet = sortrows(softTaskSet, 1);

%% 
tic
taskSetTemp = [hardTaskSet; softTaskSet];
lambdaTemp = [lambdaHard, lambdaSoft];

taskCount = hardTaskCount + softTaskCount;
A = lambdaTemp.* (taskSetTemp(:,3)./taskSetTemp(:,1))';
MaxIterations = 10000;
p = 0.1;
population = 80;
[x,fval] = ga1(hardTaskCount,softTaskCount,processorCount,population,A,MaxIterations,p)
toc
%% 

for i = 1:processorCount
    ecuTaskSet{i} = [];
end
virtprocessor=[];    %

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

%% 
hold on
grid on
for i = 1:size(pj,1)
    plot([50:10:110],pj(i,:));hold on;
end



