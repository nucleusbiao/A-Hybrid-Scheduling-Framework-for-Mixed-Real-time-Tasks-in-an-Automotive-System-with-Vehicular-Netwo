function [group] = initgroup(population,hardtaskCount,softtaskCount,processorCount)

group=zeros(hardtaskCount+softtaskCount,population);
for j = 1:population
    for i = 1:hardtaskCount
        group(i,j) = processorCount/10*rand();
    end
    for i = hardtaskCount+1:hardtaskCount+softtaskCount
        group(i,j) = (3*processorCount)/10*rand;
    end
end
    
end

