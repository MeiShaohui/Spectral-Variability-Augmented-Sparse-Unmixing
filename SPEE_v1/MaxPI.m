function PI = MaxPI(kernel)
%%%%%%%%%%%%%%%%%%calculate the Maximum distance purity index%%%%%%%%%%%%%%%
[b,num] = size(kernel);
ck = mean(kernel')';
for k = 1:num
    temp = kernel(:,k);
    sam(k) = sum(temp.*ck)/sqrt(sum(temp.^2))/sqrt(sum(ck.^2));
end          

PI = min(sam);