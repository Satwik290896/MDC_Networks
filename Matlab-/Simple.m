clear all; 
lambda=1000;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);

%very rarely plotted
med_proc = rand(90,2)

D = pproc/med_proc;