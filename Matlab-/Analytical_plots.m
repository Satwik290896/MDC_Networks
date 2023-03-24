clc;
clear all;


SINR_AWGN = [-1]; % in dB
var = 10.^(SINR_AWGN/10);%Conversion from dB


R_c = 0.4
%TT=0.1
freq = 1.5*10^9
wave_l = (3*10^8)/freq
noise = 0.1
Thresh = linspace(-30,15,300)
TT = 10.^(Thresh/10);

Mul = wave_l*3.14*((R_c)^2);

%%%  USING EQN IN REF PAPER 3 
% Get D from Simulation.
c = 1;
for m=1
   for i=1:300
    func = @(rxxx) exp((-Mul*U2(TT(i),D,rxxx,2,2,4))-((TT(i)*var(m)+TT(i)*c)*(R_c^2)*rxxx));
    func1 = @(rxxx) exp((-Mul*U2(TT(i),D,rxxx,2,2,4))-((TT(i)*var(m)+TT(i)*c)*(R_c^4)*rxxx^2));
    Int_func1(i) = integral(func1,1,Inf)
    Int_func(i) = integral(func,0,1)

    Sim_AvgPr(i) = Mul * (Int_func(i)+Int_func1(i))

   end

 figure(2)
 hold on


plot(Thresh,Sim_AvgPr,'DisplayName',"SINR:"+ num2str(SINR_AWGN(m))+"dB")

title('SINR Coverage Probability','FontSize',12)
xlabel('Threshold in dB','FontSize',10)
ylabel('P(SINR>T)','FontSize',10)
legend()
legend('Location','northeast')



if (m==1)
    writematrix(Sim_AvgPr,"SINR_185dB.txt")
    writematrix(Thresh, "Thresh.txt")
elseif (m==2)
    writematrix(Sim_AvgPr,"SINR_3dB.txt")
elseif (m==3)
    writematrix(Sim_AvgPr, "SINR_7dB.txt")
end

end


%function [out] = I1(Tt,r,alp_0,d,alp_1)
    
    % Tt Threshold > vector 
    % r : integrand 
    % r is the integrand, everything related to it should have dot operator
    % alp_0, d, alp_1 : scalars

    %out = hyp_geo((alp_0./d),1./(Tt.*(r.^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt.*(r.^(alp_0/d)))) - (r.*hyp_geo((alp_0/d),1./Tt)) + r - 1;
    
    %out = hyp_geo((alp_0/d),1/(Tt*(r^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt*(r^(alp_0/d)))) - (r*hyp_geo((alp_0/d),1/Tt)) + r - 1;

%    out = hyp_geo((alp_0/d),1./(Tt*(r.^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt*(r.^(alp_0/d)))) - (r.*hyp_geo((alp_0/d),1/Tt)) + r - 1;


%end