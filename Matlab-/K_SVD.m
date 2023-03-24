clear all; 
lambda=1000;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);

un_Array = rand(npoints,npoints);
%un_Array(un_Array<0.7) = 0;
N_const = sqrt(sum(un_Array.*un_Array,1));
N_Array = un_Array./(repmat(N_const,npoints,1));

%X1 = wmpalg('OMP',pproc(:,1),N_Array);

%X2 = wmpalg('OMP',pproc(:,2),N_Array);

%X = [X1 X2];

opt = struct('ReguParaType','SLEP','ReguAlpha',0.05);
mA = SparseCodingwithBasis(N_Array,pproc,opt);
X = cell2mat(mA)';

figure(1)
plot(pproc(:,1),pproc(:,2),'.b')

figure(2)
plot(X(:,1),X(:,2),'.b')

tmp_3D = zeros(npoints,2,npoints);

for k =1:npoints
    tmp_3D(:,:,k) = tmp_3D(:,:,k) + N_Array(:,k) * X(k,:);
end

t_s = sum(tmp_3D,3);
fx = pproc - t_s;

E_k3D = repmat(fx,1,1,npoints)+tmp_3D;

ome_3D = zeros(2,2,npoints);
[r,c] = find(X == 0);

D_new = zeros(npoints,npoints);
X_new = zeros(npoints,2);
for i = 1:npoints
    a = r == i;
    switch sum(a)
        case 1
            [m,id] = max(a);
            X_k = nonzeros(X(i,:));
            if(c(id) == 1)
                ome_3D(:,1,i) = [0;1];
                E_k = E_k3D(:,:,i)*ome_3D(:,1,i);
                [U,S,V] = svd(E_k);
                X_k = V(:,1);
                X_new(i,:) = [0 X_k(1,1)];
                D_new(:,i) = U(:,1);
            else
                ome_3D(:,1,i) = [1;0];
                 E_k = E_k3D(:,:,i)*ome_3D(:,1,i);
                [U,S,V] = svd(E_k);
                X_k = V(:,1);
                X_new(i,:) = [X_k(1,1) 0];
                D_new(:,i) = U(:,1);
            end
                            
        case 0
            ome_3D(:,:,i) = eye(2);
            [U,S,V] = svd(E_k3D(:,:,i));
            D_new(:,i) = U(:,1);
            X_new(i,:) = V(:,1);
            
        case 2
            D_new(:,i) = zeros(npoints,1);
            X_new(i,:) = zeros(1,2);
    end
end


X = X_new;
N_const = sqrt(sum(D_new.*D_new,1));
N_Array = D_new./(repmat(N_const,npoints,1));

figure(3)
plot(X_new(:,1),X_new(:,2),'.');


            