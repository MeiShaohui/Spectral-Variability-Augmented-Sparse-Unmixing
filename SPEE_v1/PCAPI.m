function PI = PCAPI(kernel)
%%%%%%%%%%%%%%%%%%%%%%%calculate the PCA purity index%%%%%%%%%%%%%%%
[COEFF, SCORE, LATENT] = PRINCOMP(kernel);
suml = sum(LATENT);
PI = LATENT(1)/suml; 