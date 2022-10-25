clear all
load jasperRidge2_R198.mat
load jasperendlist.mat

image = reshape(Y',nRow,nCol,nBand);
p = 4; 

%%%%%%%%%%%%%%%%%%%%%%%%SPEE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endlist = SPEE(image,50,5,0.86); %%5-SVD
[ReconTrain,VariaTrain] = pca_forImage(endlist,p);
library_image = ReconTrain;
library_varian = VariaTrain;

%%%%%%%%%%%%%%%%%%%%%%%%SVASU%%%%%%%%%%%%%%%%%%%%%%%%%%%
lamda1 = 5;
lamda2 = 1e4;
lamda3 = 1e4;
maxiter = 50;
[A,B,Cost] = SVASU(Y,p,maxiter,library_image,library_varian,1,lamda1,lamda2,lamda3);
A = A./repmat(sum(A),size(endlist,2),1);
Recon_image = library_image*A+library_varian*B;   

%%%%%%%%%%%%%%%%%%%%%%%%%Evaluation%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allsreimage = SRE(Y,Recon_image)
rmseimage = RMSE(Y,Recon_image);
allrmse_image = mean(mean(rmseimage))
