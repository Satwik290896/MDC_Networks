function [D] = Correlating_function(npoints,med_proc,pp,XY_disp,n1)

pproc = pp - repmat(XY_disp,npoints,1)
%figure(1)
%plot(pproc(:,1),pproc(:,2),'.');

A = zeros(n1,n1);

for i=1:n1-1
    A(i,i+1) = 1;
end

C = med_proc';
Dt = pproc';
b1 = zeros(n1,1)+0.00001;
b1(n1,1) = 0 ;
b2_t = zeros(n1,1);
b3_t = zeros(n1,1)+0.00001;
D = [];
Array = zeros(n1-2,n1);
C1 = [C;Array]; 

for i = 1:npoints
    j = rem(i,n1);
    if j==0
        j = n1;
    end
    A1 = A;
    b2 = b2_t;
    b3 = b3_t;
    b3(j) = 10;
    b2(j) = -10;

    if j~=1 
        A1(j-1,:) = [1 zeros(1,n1-1)];
    end
    disp(i);
    Xt = lsqlin(C1,[Dt(:,i);zeros(n1-2,1)],A1,b1,[],[],b2,b3);
    D = [D;Xt'];
end

%D = pproc/med_proc;

%summation_of_distances
