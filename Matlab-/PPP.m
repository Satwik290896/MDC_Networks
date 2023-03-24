clear all; 
lambda=3000;
figure
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);
p1 = rand;
p2 = rand;
receiver = [p1 p2];
plot(pproc(:,1)',pproc(:,2),'.g',p1,p2,'.r')
hold on
r = 0.1;
circle(p1,p2,r);

e_d = pproc-repmat(receiver,npoints,1);

e = e_d.*e_d;
d = sqrt(sum(e,2));

b = d(d<0.05);
pi = randi([1,length(b)]);
where = d == b(pi,1);
id = find(where);
emitter = [pproc(id,1) pproc(id,2)];
plot(pproc(id,1),pproc(id,2),'.b');

dist_e_r = sqrt(sum((emitter-receiver).^2));
beta = 10^-30;
P_e = exp(-npoints*(4/(2*sqrt(2)))*10*(dist_e_r^3)*beta);







 