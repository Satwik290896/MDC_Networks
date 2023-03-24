
med_temp = sum((med_proc.^2),2);
S = (D.^2)*med_temp;
for i = 1:npoints
    for j = 1:npoints
        sum = 0;
        sum = sum + (D(i,j)^2)*(med_proc(j,1)^2+med_proc(j,2)^2);
    end
    S = [S;sqrt(sum)];
end
P = 0;
for i = 1:npoints
    for j = i+1:npoints
        P = P + 2*D(1,i)*D(1,j)*(med_proc(i,1)*med_proc(j,1)+med_proc(i,2)*med_proc(j,2));
    end
end
    
P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;
p11 = 0.4+rand(100,1)*0.2;
p22 = 0.4+rand(100,1)*0.2;
receiver = [p11 p22];