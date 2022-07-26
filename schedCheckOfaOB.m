function [schedFlag, deadlineLOArray] = schedCheckOfaOB(taskSet, initialOB)
%taskset(i) = [period, deadline, wcetLO, wcetHI, criticalityLevel];

serviceTask = rtcfsu(1);
serviceTaskLO = rtcaffine(serviceTask, 1,  initialOB);
serviceTaskHI = serviceTask;
schedFlag = 1;
approxLength = 3*max(taskSet(:,1));
deadlineLOArray = [];
% figure
% hold on
% rtcplot(serviceTaskLO, approxLength, 'k')
% rtcplot(serviceTaskHI, 'k', 'LineWidth',2)
for i = 1:size(taskSet, 1)
    [deadlineLO, schedFlag] = FindMinimalVirtualDeadline(taskSet(i,:), serviceTaskLO, serviceTaskHI);
    deadlineLOArray = [deadlineLOArray, deadlineLO];
    if schedFlag == 0
        break
    else
        wbfLOi = rtctimes(rtcpjdu(taskSet(i,1), 0, 0), taskSet(i,3));
        wbfLOi = rtcapproxs(wbfLOi, approxLength, 0, 1);
%         rtcplot(wbfLOi, 'r')
        serviceTaskLO = rtcmaxconv(rtcminus(serviceTaskLO, wbfLOi), 0);
%         rtcplot(serviceTaskLO, approxLength, 'k')
        if taskSet(i, 5) == 1
            baseHIi = rtctimes(rtcpjdu(taskSet(i,1), 0, 0), taskSet(i,4));
            wbfHIi = rtcaffine(baseHIi, 1, -deadlineLO);
            wbfHIi = rtcapproxs(wbfHIi, approxLength, 0, 1);
%             rtcplot(wbfHIi, 'r', 'LineWidth',2)
            serviceTaskHI = rtcmaxconv(rtcminus(serviceTaskHI, wbfHIi), 0);
%             rtcplot(serviceTaskHI, 'k', 'LineWidth',2)
        end
    end
end
%pause