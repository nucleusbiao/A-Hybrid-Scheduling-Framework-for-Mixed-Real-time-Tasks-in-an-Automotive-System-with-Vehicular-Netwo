function fit = fitness(group,A,processorCount)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
fit=zeros(1,size(group,2));

for j = 1:size(group,2)
    n = zeros(size(A,1),1);
    u = zeros(size(A,1),1);
    for i = 1:size(group,1)
        if group(i,j) < processorCount/10
            fit(j) = fit(j) + 1;
            n(ceil(group(i,j)*10)) = n(ceil(group(i,j)*10)) + 1;
            u(ceil(group(i,j)*10)) = u(ceil(group(i,j)*10)) + A(ceil(group(i,j)*10),i);
        end     
    end
    for i = 1:length(u)
        if u(i) > n(i)*(2^(1/n(i)) - 1)
            fit(j) = 0;
            break
        end
    end
end
end

