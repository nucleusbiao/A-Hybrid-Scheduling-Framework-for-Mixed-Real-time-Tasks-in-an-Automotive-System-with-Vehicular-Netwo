function [schFlag,processorID,taskSet,pos] = schTestOnline(taskSet,insertTask,lambdaSoftLeft, t)

processorID = 0;
pos = 0;
ecuUtilization = zeros(size(taskSet,1),1);
for i = 1:size(taskSet,1)
    for j = 1:size(taskSet,2)
        if ~isempty(taskSet{i,j})
        ecuUtilization(i) = ecuUtilization(i) + taskSet{i,j}{5} * taskSet{i,j}{2}/taskSet{i,j}{1};
        end
    end
end
[U,processorID] = min(ecuUtilization);
taskSetTemp = taskSet(processorID,:);
for i = 1:size(taskSet,2)
    if isempty(taskSetTemp{i})
        taskNum = i - 1;
        break
    end
end
approxLength = 3*insertTask{1};
for m = 1:size(taskSetTemp,2)
    if ~isempty(taskSetTemp{m})
        if 3*taskSetTemp{m}{1} > approxLength
            approxLength = 3*taskSetTemp{m}{1};
        end
    end
end

for pos = 0:taskNum
    serviceECU = rtcfsu(1);
    if pos == 0
        taskSetTemp{taskNum+1} = insertTask;
    else
        taskSetTemp{taskNum-pos+1} = insertTask;
        taskSetTemp(taskNum-pos+2:taskNum+1) = taskSet(processorID,taskNum-pos+1:taskNum);
    end
    for i = 1:taskNum+1
        schFlag = true;
        if ~isempty(taskSetTemp{i})
            if i == taskNum - pos + 1
                wbfTask = rtctimes(rtcpjdu(taskSetTemp{i}{1}, 0, 0),taskSetTemp{i}{2}*lambdaSoftLeft(insertTask{10},processorID));
            else                
                wbfTask = rtctimes(rtcpjdu(taskSetTemp{i}{1}, 0, 0),taskSetTemp{i}{2}*taskSetTemp{i}{5});
            end
            if t-taskSetTemp{i}{7} ~= 0
                wbfTask = rtcaffine(wbfTask, 1,  -(t-taskSetTemp{i}{7}));
            end
            if taskSetTemp{i}{8} ~=0
                c = rtccurve([[0 -taskSetTemp{i}{3} 0];[1 -taskSetTemp{i}{3} 0]]);
            else
                c = rtccurve([[0 -taskSetTemp{i}{2} 0];[1 -taskSetTemp{i}{2} 0]]);
            end
            wbfTask = rtcplus(wbfTask, c);
            wbfTask = rtcapproxs(wbfTask, approxLength, 0, 1);  % + (t-taskSetTemp{i}{7})
            demandTask = rtcaffine(wbfTask, 1,  taskSetTemp{i}{4});
            if rtceq(serviceECU, rtcmax(serviceECU, demandTask))
                serviceECU = rtcmaxconv(rtcminus(serviceECU, wbfTask), 0);
            else
                schFlag = false;
                break
            end
        end
    end
    if schFlag == true
        break
    end
end

if schFlag == true
    pos = taskNum - pos + 1;
end


end











% for i = 1:size(taskSet,1)
%     serviceECU = [serviceECU rtcfsu(1)];
% end
% for i = 1:size(taskSet,1)
%     hold off;
%     approxLength = 0;
%     for m = 1:size(taskSet,2)
%         if ~isempty(taskSet{i,m})
%             if 3*taskSet{i,m}{1} > approxLength
%                 approxLength = 3*taskSet{i,m}{1};
%             end
%         end
%     end
%     for j = 1:size(taskSet,2)
%         hold off
%         if ~isempty(taskSet{i,j})           
%             wbfTask = rtctimes(rtcpjdu(taskSet{i,j}{1}, 0, 0), taskSet{i,j}{2}*taskSet{i,j}{5});
%             if t-taskSet{i,j}{7} ~= 0
%                 wbfTask = rtcaffine(wbfTask, 1,  -(t-taskSet{i,j}{7}));
%             end
%             if taskSet{i,j}{8} ~=0
%                 c = rtccurve([[0 -taskSet{i,j}{3} 0];[1 -taskSet{i,j}{3} 0]]);
%             else
%                 c = rtccurve([[0 -taskSet{i,j}{2}*taskSet{i,j}{5} 0];[1 -taskSet{i,j}{2}*taskSet{i,j}{5} 0]]);
%             end
%             wbfTask = rtcplus(wbfTask, c);
%             wbfTask = rtcapproxs(wbfTask, approxLength+(t-taskSet{i,j}{7}) , 0, 1);  %+ (t-taskSet{i,j}{7});
%             rtcplot(wbfTask, 'b--',10000);hold on;
%             demandTask = rtcaffine(wbfTask, 1,  taskSet{i,j}{4});
%             rtcplot(demandTask, 'g--',10000);hold on;
%             rtcplot(serviceECU(i), 'r--',10000);hold on;
% 
%             if rtceq(serviceECU(i), rtcmax(serviceECU(i), demandTask))
%                 serviceECU(i) = rtcmaxconv(rtcminus(serviceECU(i), wbfTask), 0);
% 
%             else
%                 
%                 display('error')
%                 schFlag = false;
% 
%             end
%         end
%     end
% end
% maxserviceECU = serviceECU(1);
% d=100;
% hold off
% for i = 1:size(taskSet,1)
%     rtcplot(serviceECU(i), 'r--',10000);hold on;
%     maxserviceECU = rtcmax(maxserviceECU,serviceECU(i));
% end
% for i = 1:size(taskSet,1)
%     dist = rtcplotv(maxserviceECU, serviceECU(i));
%     if dist < d
%         d = dist;
%         processorID = i;
%     end
% end
