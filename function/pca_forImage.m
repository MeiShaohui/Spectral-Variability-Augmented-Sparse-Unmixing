function [ReconTrain,VariaTrain]=pca_forImage(DataTrain,k)
%% function: 求解某矩阵的PCA估计值
%% input
%DataTrain: 待映射样本 bands*samples
%T：阈值   
%k：个数

%% output
%ReconTrain：   Data_train提纯光谱                
%VariaTrain：   Data_train光谱变化                

[nb samples] = size(DataTrain);
%归一化
varV = std(DataTrain,0,2);
Meanvalue = mean(DataTrain');
MeanData = repmat(Meanvalue',[1 samples]);
DeMeanData = DataTrain-MeanData;
CovData = DeMeanData*DeMeanData'/samples;
[U,D]=eig(CovData);                    %协方差矩阵的特征值D,特征向量U
[C,I]=sort(diag(D),'descend');         %特征值从大到小排序

%阈值选取前k个
% Allsum = sum(C);
% sumC = 0;
% for k = 1:size(C,1)
%     sumC = sumC+C(k);
%     if sumC/Allsum > T
%         break;
%     end
% end
Uk = U(:,I(1:k));
Uend = U(:,I(k+1:end));

Coeff_k = Uk'*DeMeanData;                  
ReconData = Uk*Coeff_k;
ReconTrain = ReconData+MeanData;
Coeff_end = Uend'*DeMeanData;
VariaTrain = Uend*Coeff_end+MeanData;

end