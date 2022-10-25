function [value] = SRE(matrix_real,matrix_est)
%% SRE signal-to-reconstruction error
%SRE(dB)=10log10(E[||x||_2 ^2]/E[||x-x^||_2 ^2])
% x-- true abundance
% x^--estimated abundance
% larger is better
% flag = 0-alpha  1-pixel reconstruct 
[nr,nl,nb] = size(matrix_est);

if nb == 1
	for j = 1:nl
		X_est = matrix_est(:,j);
		X_real = matrix_real(:,j);
		Eu(j) = norm(X_real,2)^2;
		Ed(j) = norm(X_real-X_est,2)^2;
	end
else
	for i = 1:nb
		X_est = reshape(matrix_est(:,:,i),nr*nl,1);
		X_real = reshape(matrix_real(:,:,i),nr*nl,1);
		Eu(i) = norm(X_real,2)^2;
		Ed(i) = norm(X_real-X_est,2)^2;
    end
end
value = 10*log10(mean(Eu)/mean(Ed));

end