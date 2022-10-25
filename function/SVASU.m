function [A,B,Cost] = VCSUV2(R,numofend,maxiter,W,V,lamda0,lamda1,lamda2,lamda3)
%variation concerned sparse unmixing
% min ||R-MA||_F + lamda1*||R-MA-VB||_F + lamda2*||A||_{2,1} + lamda3*||B||_F

[nb samples] = size(R);
%初始化
rand('state',0)
A = rand(size(W,2),samples);
rand('state',1)
B = rand(size(V,2),samples);

%和为一
dita = 20;%mean(mean(R));
tR = [R; dita*ones(1,samples)];
tW = [W; dita*ones(1,size(W,2))];
tV = [V; dita*ones(1,size(V,2))];%20
tic

for iter=1:maxiter,

    %% update A
     for k = 1 : size(A,1)
        D(k,k) = 1/norm(A(k,:),2);
     end
     multiA = (lamda0+lamda1)*tW'*tR-lamda1*tW'*tV*B;
     divisA = (lamda0+lamda1)*tW'*tW*A+lamda2*D*A;
     A = A.*multiA./divisA;
     A=A./repmat(sum(A),size(A,1),1);
	
    %% update B  
     multiB = lamda1*tV'*tR-lamda1*tV'*tW*A;
     divisB = lamda1*tV'*tV*B+lamda3*B;
     B = B.*multiB./divisB;
      
    %%  cost   
    term0(iter) = norm((R-W*A),'fro');
    term1(iter) = norm((R-W*A-V*B),'fro');
    term2(iter) = L21(A);
    term3(iter) = norm(B,'fro');

    Cost(iter) = lamda0*term0(iter) + lamda1*term1(iter) + lamda2*term2(iter) + lamda3*term3(iter);
    %fprintf('iter = %d, loss = %d\n',iter,Cost(iter));
    
end
end