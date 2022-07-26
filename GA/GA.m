clear
periodArray = [20 25 40 50 80 100 200 250 400 500 800 1000];
load('F:\work\ÂÛÎÄ-1\1\processor20task200.mat')
processorCount = 20;
hardTaskCount = 100;
softTaskCount = 100;
% vectU = UUniFast(40,9);
% for i = 1:hardTaskCount
%     period = periodArray(randi(12));                                    
%     hardTaskSet(i,:) = [period, period, ceil(vectU(randi(40))*period)];      %å‘å¸ƒå‘¨æœŸ  æœ?æœŸé™  
% end
% 
% for i = 1:softTaskCount
%     period = periodArray(randi(12));
%     softTaskSet(i,:) = [period, period, ceil(vectU(randi(40))*period)];
% end
% 
% lambdaHard = (5+randi(10, processorCount, hardTaskCount))/10;
% lambdaSoft = (5+randi(10, processorCount, softTaskCount))/10;
% 
% hardTaskSet = sortrows(hardTaskSet, 1);
% softTaskSet = sortrows(softTaskSet, 1);
GAn = zeros(20,1);
GAt = zeros(20,1);
for k = 1:20
    %     hardTaskSet0 = eval(['hardTaskSet',num2str(k)]);
    %     softTaskSet0 = eval(['softTaskSet',num2str(k)]);
    %     lambdaHard0 = eval(['lambdaHard',num2str(k)]);
    %     lambdaSoft0 = eval(['lambdaSoft',num2str(k)]);
    hardTaskSet0 = hardTaskSet{k};
    softTaskSet0 = softTaskSet{k};
    lambdaHard0 = lambdaHard{k}(1:processorCount,:);
    lambdaSoft0 = lambdaSoft{k}(1:processorCount,:);
    
    hardTaskSet0 = sortrows(hardTaskSet0, 1);
    softTaskSet0 = sortrows(softTaskSet0, 1);
%% 

taskSetTemp = [hardTaskSet0; softTaskSet0];
lambdaTemp = [lambdaHard0, lambdaSoft0];

taskCount = hardTaskCount + softTaskCount;
A = [];
for i = 1:size(lambdaTemp,1)
    A(i,:) = lambdaTemp(i,:).* (taskSetTemp(:,3)./taskSetTemp(:,1))';
end
MaxIterations = 50000;
p = 0.1;
population = 100;
t1 = clock;
[x,fval] = ga1(hardTaskCount,softTaskCount,processorCount,population,A,MaxIterations,p);
t2 = clock;
%% 
GAn(k) = fval;
GAt(k) = etime(t2,t1);
for i = 1:processorCount
    ecuTaskSet{i} = [];
end
virtprocessor=[];    %è™šæ‹Ÿå¤„ç†å™?

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
%% 




