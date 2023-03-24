clear all;
%Formulae g(t)%
emp = [];
Thre = linspace(-30,120,500);

for xxxx = Thre
    emp = [emp (1-(1/(1+exp(-xxxx))))];
end

figure(1)
plot(Thre,emp);

fileID = fopen('Values_T_TrueAWGN_-5dB.txt.txt','a+');
fprintf(fileID,'%f\n',emp);
fclose(fileID);