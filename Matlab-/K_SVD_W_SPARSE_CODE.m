clear all; 
lambda=100;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);

figure(1)
plot(pproc(:,1),pproc(:,2),'.b');
un_Array = rand(npoints,npoints);
un_Array(un_Array<0.9) = 0;
N_const = sqrt(sum(un_Array.*un_Array,1));
N_Array = un_Array./(repmat(N_const,npoints,1));

X = DL_KSVD_test(N_Array, pproc);

figure(2)
plot(X(:,1),X(:,2),'.b');