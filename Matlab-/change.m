clear all; 
lambda=3000;
figure
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);


plot(pproc(:,1)',pproc(:,2),'.g',p1,p2,'.r')
hold on

r = 0.1;
circle(p1,p2,r);

b = d(d<0.03);
pi = randi([1,length(b)]);
where = d == b(pi,1);
id = find(where);

plot(pproc(id,1),pproc(id,2),'.b');

d_er =  sqrt((pproc(id,1)-p1)^2+(pproc(id,2)-p2)^2)