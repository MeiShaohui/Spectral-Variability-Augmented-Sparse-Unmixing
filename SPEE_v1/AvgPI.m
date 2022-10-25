function PI = AvgPI(kernel)
%%%%%%%%%%%%%%%%%%calculate the average distance purity index%%%%%%%%%%%%%%%
[b,num] = size(kernel);
ck = mean(kernel')';
for k = 1:num
    temp = kernel(:,k);
    sam(k) = sum(temp.*ck)/sqrt(sum(temp.^2))/sqrt(sum(ck.^2));
end          

PI = mean(sam);