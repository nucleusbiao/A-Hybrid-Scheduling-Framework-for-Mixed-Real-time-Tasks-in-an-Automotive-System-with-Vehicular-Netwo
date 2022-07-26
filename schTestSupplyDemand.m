function schFlag = schTestSupplyDemand(taskArray)


serviceECU = rtcfsu(1);
approxLength = 3*max(taskArray(:,1));

schFlag = true;

for i = 1:size(taskArray, 1)
    hold off
    wbfTask = rtctimes(rtcpjdu(taskArray(i,1), 0, 0), taskArray(i,3)*taskArray(i,4));
    wbfTask = rtcapproxs(wbfTask, approxLength, 0, 1);
%     rtcplot(wbfTask, 'b--',10000);hold on;
    demandTask = rtcaffine(wbfTask, 1,  taskArray(i,2));
%     rtcplot(demandTask, 'g--',10000);hold on;
%     rtcplot(serviceECU, 'r--',10000);hold on;
    if rtceq(serviceECU, rtcmax(serviceECU, demandTask))
        serviceECU = rtcmaxconv(rtcminus(serviceECU, wbfTask), 0);
    else
        schFlag = false;
        break
    end
end