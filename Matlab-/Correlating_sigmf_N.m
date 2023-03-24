clear all;
%Formulae g(t)%
emp = [];
Thre = linspace(-30,120,500);

for xxxx = Thre
    emp = [emp (1-sig(exp(xxxx/10.0)))];
end

figure(1)
plot(Thre,emp);

fileID = fopen('Values_T_Sigmoid_N.txt.txt','a+');
fprintf(fileID,'%f\n',emp);
fclose(fileID);

function [h] = sig(t)
    h = 1/(1+exp(-t))
end