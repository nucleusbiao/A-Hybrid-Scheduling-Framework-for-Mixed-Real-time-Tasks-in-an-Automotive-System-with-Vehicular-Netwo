function [single] = mutation(single,processorCount)

    pos = ceil(rand()*size(single,1));
    single(pos) = (processorCount + 5)/10*rand;

end

