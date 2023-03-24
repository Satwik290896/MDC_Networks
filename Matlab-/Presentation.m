clear all;
lambda=200;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2)-0.5;
XY_dar = 0.5*rand(1,2);
pp = pproc - repmat(XY_dar,npoints,1);
n1 = 2;
med_proc = rand(8,2);

figure(1)
plot(med_proc(:,1),med_proc(:,2),'.');
hold on;
plot(0,0,'.r')

figure(2)
plot(pp(:,1),pp(:,2),'.');
hold on;
plot(0,0,'.r');