function PI = MEI2(pixel)
%%%%%%%%%%%%%%%%%%%%%%%calculate the MEI purity index%%%%%%%%%%%%%%%
%% to every pixel of the kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%caculate the cumulative distance to every pixel
[b,num] = size(pixel);
sam = zeros(1,num);
for k1=1:num
    ck = pixel(:,k1);
    for k2=1:num
        if k1~=k2
            temp = pixel(:,k2);
            sam(k1) = sam(k1)+ sum(temp.*ck)/sqrt(sum(temp.^2))/sqrt(sum(ck.^2));
        end
    end
end

[tt,locmax] = max(sam);
dilation = pixel(:,locmax);
[tt,locmin] = min(sam);
erosion = pixel(:,locmin);      
% set the MEI to the max pixel in the kernel
PI = sum(dilation.*erosion)/sqrt(sum(dilation.^2))/sqrt(sum(erosion.^2));