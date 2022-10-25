function rmsevalue = RMSE(a,b)
%%%%%%%%%%%%RMSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% caculate the root mean square error (RMSE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2008-01-10

[row,line,c] = size(a);

temp1 = zeros(1,c);
temp2 = zeros(1,c);
for i = 1:row
    for j = 1:line
        temp1(1:c) = a(i,j,1:c);
        temp2(1:c) = b(i,j,1:c);
        rmsevalue(i,j) = sum((temp1-temp2).^2);
    end
end

rmsevalue = sqrt(rmsevalue./row./line);
end