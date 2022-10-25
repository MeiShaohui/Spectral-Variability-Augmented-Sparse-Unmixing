function PI = SVDPI(kernel)
%%%%%%%%%%%%%%%%%%%%%%%%%calculate the SVD purity index%%%%%%%%%%%%%%%
[s,v,d] = svd(kernel);
sumv = 0;
m = rank(v);
for k=1:m
    sumv = sumv+v(k,k);
end
PI = v(1,1)/sumv;