function PI = MEI1(kernel)
%%%%%%%%%%%%%%%%%%%%%%%calculate the MEI purity index%%%%%%%%%%%%%%%
%% to the average centre of the kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b,num] = size(kernel);
ck = mean(kernel')';
for k = 1:num
    temp = kernel(:,k);
    sam(k) = sum(temp.*ck)/sqrt(sum(temp.^2))/sqrt(sum(ck.^2));
end          
[tt,locmax] = max(sam);
dilation = kernel(:,locmax);
[tt,locmin] = min(sam);
erosion = kernel(:,locmin);
PI = sum(dilation.*erosion)/sqrt(sum(dilation.^2))/sqrt(sum(erosion.^2));