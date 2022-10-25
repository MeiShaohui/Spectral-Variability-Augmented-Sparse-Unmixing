function endlist = SPEE(pixel,numofendmember,kind,gama)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function: spatial purity based endmember extraction algorithm
%% endlist: the final endmember
%% pixel: the initial image
%% numofendmember: the num of endmember needed to be generate 
%% kernel: 3x3 window
%% kind: 1-average distance 2-maximum distance 3-AMEE(to average) 4-AMEE(to
%%          every pixel) 5-SVD 6-PCA
%% gama: threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row,line,b] = size(pixel);
pixel1 = pixel;

image_index = zeros(row,line);
for i = 2:row-1
    for j = 2:line-1
        kernel(1:b,1) = pixel(i-1,j-1,1:b);
        kernel(1:b,2) = pixel(i-1,j,1:b);
        kernel(1:b,3) = pixel(i-1,j+1,1:b);
        kernel(1:b,4) = pixel(i,j-1,1:b);
        kernel(1:b,5) = pixel(i,j,1:b);
        kernel(1:b,6) = pixel(i,j+1,1:b);
        kernel(1:b,7) = pixel(i+1,j-1,1:b);
        kernel(1:b,8) = pixel(i+1,j,1:b);
        kernel(1:b,9) = pixel(i+1,j+1,1:b);
                
        switch kind
            case 1
                image_index(i,j) = AvgPI(kernel);
            case 2
                image_index(i,j) = MaxPI(kernel);
            case 3
                image_index(i,j) = MEI1(kernel);
            case 4
                image_index(i,j) = MEI2(kernel);
            case 5
                image_index(i,j) = SVDPI(kernel);
            case 6
                image_index(i,j) = PCAPI(kernel);
            otherwise
        end
            
        if image_index(i,j)>gama         % use the average pixel of pixels in a kernel
            temp = mean(kernel');
            pixel1(i,j,1:b) = temp(1:b);
        end
    end
end

for i = 1:row
    for j = 1:line
        if image_index(i,j)>gama
            image_index1(i,j) = 1;
        else
            image_index1(i,j) = 0;
        end
    end
end

%%%%%%%%% show the SN purity map to determine suitable threshold
figure
imshow (image_index1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%spatial combination%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dita = 2;
index_end = 1;
flag_check = zeros(row,line);
for i = 1:row
    for j = 1:line
        if ((flag_check(i,j)==0) & (image_index1(i,j)==1))
            flag_check(i,j)=1;
            list = zeros(2,1);
            num = 1;        %队列的元素数目
            list(1,num) = i;
            list(2,num) = j;
            index = 1;      %队列的头元素下标

            while index<=num
                s = list(1,index);
                t = list(2,index);
                index = index+1;
                for k=-1:1:1
                    for l = -1:1:1
                        if((k~=0)||(l~=0))
                            if (image_index1(s+k,t+l)==1)&(flag_check(s+k,t+l)==0)
                                num = num+1;
                                list(1,num) = s+k;
                                list(2,num) = t+l;
                                flag_check(s+k,t+l)=1;
                            end
                        end
                    end
                end
            end
            if num==1
                endlist(1:b,index_end) = pixel1(i,j,1:b);          %光谱均值
                coordinate(1,index_end) = i;
                coordinate(2,index_end) = j;
                index_end = index_end+1;
            else
                %%%%%%%%%%%%%%%%%%%%%%求出连通域内所有像元的均值，作为纯净像元%%%%%%%%%%%%%%%%%%%%%%%%%
                sumend = zeros(b,1);
                for p = 1:num
                    for q = 1:b
                        sumend(q) = sumend(q) + pixel1(list(1,p),list(2,p),q);
                    end
                end
                endlist(:,index_end) = sumend/num;          %光谱均值
                avcor = mean(list');                        %坐标均值
                coordinate(1,index_end) = avcor(1);
                coordinate(2,index_end) = avcor(2);
                index_end = index_end+1;
            end
        end
    end
end
clear pixel1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show the location of all spatially independent candidates 
figure
temp = pixel(:,:,2:2:6);
imshow(uint8(temp));
c = index_end-1;
for i=1:c
    hold on;
    plot(coordinate(2,i),coordinate(1,i),'+');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  combination by the spectral relation %%%%%%%%
% every time combine two most similar spectrum
flag = 1;
sam = zeros(c,c);
while c>numofendmember    %
    if flag==1
        for i=1:c
            for j=1:c
                if i~=j
                    sam(i,j) = endlist(:,i)'*endlist(:,j)/sqrt(endlist(:,i)'*endlist(:,i))/sqrt(endlist(:,j)'*endlist(:,j));
                end
            end   
        end
        flag=0;
    end
    [m,locx] = max(sam);
    [mm,locy] = max(m);
    %%%% combine the locx and locy
    ml = locx(locy);
    mr = locy;
    endlist(:,ml)=(endlist(:,ml)+endlist(:,mr))/2;
    for i=1:c
        if i~=ml
            sam(ml,i) = endlist(:,i)'*endlist(:,ml)/sqrt(endlist(:,i)'*endlist(:,i))/sqrt(endlist(:,ml)'*endlist(:,ml));
            sam(i,ml) = sam(ml,i);
        end
    end
    endmember(:,1:mr-1) = endlist(:,1:mr-1);
    endmember(:,mr:c-1) = endlist(:,mr+1:c);
    sam(:,mr:c-1) = sam(:,mr+1:c);
    sam(mr:c-1,:) = sam(mr+1:c,:);    
    c=c-1;
    sam = sam(1:c,1:c);
    endlist = endmember(:,1:c);
end
